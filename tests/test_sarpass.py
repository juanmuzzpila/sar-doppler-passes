"""Checks that do not rely on the code under test agreeing with itself."""

import math
from datetime import datetime, timezone

import numpy as np
import pytest
from sgp4.api import Satrec

from sarpass import geo, sar
from sarpass.demo import demo_constellation, sso_inclination_deg
from sarpass.pipeline import build_output
from sarpass.sar import Config, Target
from sarpass.tle import checksum, make_tle, parse_tle_text

START = datetime(2026, 11, 18, 12, 0, 0, tzinfo=timezone.utc)
VALENCIA = Target("Valencia", 39.4699, -0.3763)


# ---------------------------------------------------------------- geodesy
def test_geodetic_roundtrip():
    lat = np.array([-60.0, 0.0, 39.4699, 89.0])
    lon = np.array([-170.0, 0.0, -0.3763, 120.0])
    h = np.array([0.0, 500.0, 570e3, 35786e3])
    lat2, lon2, h2 = geo.ecef_to_geodetic(geo.geodetic_to_ecef(lat, lon, h))
    assert np.allclose(lat, lat2, atol=1e-8)
    assert np.allclose(lon, lon2, atol=1e-8)
    assert np.allclose(h, h2, atol=1e-3)


def test_known_ecef_point():
    # Equator / prime meridian at sea level is exactly (a, 0, 0)
    assert np.allclose(geo.geodetic_to_ecef(0.0, 0.0, 0.0), [geo.A, 0.0, 0.0])


# -------------------------------------------------------------------- TLE
def test_synthetic_tle_is_valid_and_has_expected_orbit():
    tle = make_tle("X", 99999, START, 570.0, 97.7, 10.0, 0.0)
    assert len(tle.line1) == len(tle.line2) == 69
    assert checksum(tle.line1) == int(tle.line1[68])
    assert checksum(tle.line2) == int(tle.line2[68])
    assert parse_tle_text(f"{tle.name}\n{tle.line1}\n{tle.line2}")[0].norad == 99999
    sat = Satrec.twoline2rv(tle.line1, tle.line2)
    _, r, _ = sat.sgp4(sat.jdsatepoch, sat.jdsatepochF)
    alt = math.sqrt(sum(c * c for c in r)) - 6378.137
    assert 540 < alt < 600
    period_min = 2 * math.pi / sat.no_kozai
    assert 95.5 < period_min < 96.5  # ~96 min at 570 km


def test_bad_checksum_is_rejected():
    tle = make_tle("X", 99999, START, 570.0, 97.7, 10.0, 0.0)
    bad = tle.line1[:68] + str((int(tle.line1[68]) + 1) % 10)
    assert parse_tle_text(f"X\n{bad}\n{tle.line2}") == []


def test_sso_inclination():
    assert 97.5 < sso_inclination_deg(570.0) < 97.9


# ---------------------------------------------------- geometry / Doppler
@pytest.fixture(scope="module")
def run():
    tle = demo_constellation(START)[1]
    t, r, v = sar.propagate_ecef(tle, START, 24.0, 1.0)
    g = sar.geometry(t, r, v, VALENCIA)
    return tle, g


def full_pass(run, min_max_elevation=0.0):
    """First pass that starts and ends inside the window (not cut by its edges)."""
    tle, g = run
    for p in sar.find_passes(g, tle, Config()):
        if p.i0 > 0 and p.i1 < len(g.t) - 1 and g.elevation[p.i0 : p.i1 + 1].max() >= min_max_elevation:
            return p
    raise AssertionError("no complete pass in 24 h")


def test_satellite_altitude_is_leo(run):
    _, g = run
    assert 540e3 < g.sat_alt.min() and g.sat_alt.max() < 600e3


def test_range_rate_matches_finite_difference_of_range(run):
    """rdot from velocity must equal d(range)/dt from positions (validates TEME->ECEF incl. Earth rotation)."""
    _, g = run
    fd = np.gradient(g.range_m, g.t)
    sel = g.elevation > 5.0
    assert sel.sum() > 100
    assert np.max(np.abs(fd[sel] - g.range_rate[sel])) < 1.0  # m/s


def test_doppler_zero_crossing_at_closest_approach(run):
    _, g = run
    p = full_pass(run)
    seg = slice(p.i0, p.i1 + 1)
    i_tca = p.i0 + int(np.argmin(g.range_m[seg]))
    # rdot changes sign within +-2 s of the minimum-range sample (1 s grid, |d rdot/dt| ~ 70 m/s^2)
    assert np.sign(g.range_rate[i_tca - 2]) != np.sign(g.range_rate[i_tca + 2])
    # sign: approaching first (f_D > 0), receding after (f_D < 0)
    f = sar.doppler_hz(g.range_rate, 9.65e9)
    assert f[p.i0] > 0 > f[p.i1]


def test_doppler_sign_convention():
    assert sar.doppler_hz(-1000.0, 9.65e9) > 0  # approaching -> positive
    assert sar.doppler_hz(+1000.0, 9.65e9) < 0


def test_doppler_at_horizon_matches_closed_form(run):
    """Near-overhead pass at elevation 0: |rdot| = v * cos(90deg - nadir_angle), sin(nadir) = Re/(Re+h)."""
    _, g = run
    p = full_pass(run, min_max_elevation=70.0)
    h = g.sat_alt[p.i0]
    # ignore Earth rotation (~0.4 km/s cross-term), so allow 10 %
    sat_speed = math.sqrt(398600.4418e9 / (geo.A + h))
    nadir = math.asin(geo.A / (geo.A + h))
    rdot_expected = sat_speed * math.cos(math.pi / 2 - nadir)
    assert abs(abs(g.range_rate[p.i0]) - rdot_expected) / rdot_expected < 0.10


def test_incidence_is_zero_when_satellite_is_overhead():
    tgt = geo.geodetic_to_ecef(39.0, -0.4, 0.0)
    up = geo.up_vector(39.0, -0.4)
    r = tgt + 570e3 * up
    v = np.array([[7000.0, 0.0, 0.0]])
    g = sar.geometry(np.array([0.0]), r[None, :], v, Target("t", 39.0, -0.4))
    assert g.incidence[0] < 1e-6
    assert abs(g.range_m[0] - 570e3) < 1e-3
    assert g.off_nadir[0] < 0.2  # ellipsoid normal vs geocentric difference is negligible


def test_incidence_window_logic(run):
    tle, g = run
    cfg = Config(inc_min_deg=15.0, inc_max_deg=50.0)
    for p in sar.find_passes(g, tle, cfg):
        for a, b in p.imaging:
            assert g.incidence[a : b + 1].min() >= 15.0 - 1e-9
            assert g.incidence[a : b + 1].max() <= 50.0 + 1e-9


# ------------------------------------------------------------- pipeline
def test_pipeline_output_schema_and_consistency():
    tles = demo_constellation(START)
    out = build_output(tles, sar.MISSIONS["iceye"], VALENCIA, Config(), START, "demo")
    assert out["source"] == "demo"
    assert out["mission"]["freq_hz"] == 9.65e9
    assert abs(out["mission"]["wavelength_cm"] - 3.107) < 0.01
    assert len(out["passes"]) > 0
    ids = [p["id"] for p in out["passes"]]
    assert ids == list(range(len(ids)))
    for p in out["passes"]:
        n = len(p["t_s"])
        assert all(len(p[k]) == n for k in ("range_km", "range_rate_ms", "incidence_deg", "sat_lat", "sat_lon"))
        assert p["imaging_windows_s"]
        assert p["side"] in ("left", "right")
