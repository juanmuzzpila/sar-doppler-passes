"""SAR pass geometry and Doppler history for a fixed ground target.

Model (deliberately simple, every step is checkable):
  * SGP4 propagates the TLE -> TEME state -> ECEF state (Earth-fixed).
  * A stationary ground target has zero velocity in ECEF, so the relative
    velocity is just the satellite's ECEF velocity.
  * Radial velocity  rdot = u . v,   u = unit vector target -> satellite.
  * Two-way (monostatic) Doppler  f_D = -2 rdot / lambda = -2 f0 rdot / c.
      approaching (rdot < 0)  ->  f_D > 0
  * Incidence angle = angle between u and the local ellipsoid normal at the
    target (flat terrain). Elevation = 90 deg - incidence.
  * Off-nadir (look) angle is measured at the satellite between the line of
    sight and the geodetic nadir.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone

import numpy as np
from sgp4.api import Satrec, jday

from . import geo
from .tle import TLE


@dataclass(frozen=True)
class Mission:
    key: str
    label: str
    celestrak_query: str
    freq_hz: float  # carrier frequency
    band: str

    @property
    def wavelength_m(self) -> float:
        return geo.C / self.freq_hz


MISSIONS = {
    "iceye": Mission("iceye", "ICEYE", "ICEYE", 9.65e9, "X"),
    "sentinel1": Mission("sentinel1", "Sentinel-1", "SENTINEL-1", 5.405e9, "C"),
}


@dataclass
class Target:
    name: str
    lat: float
    lon: float
    h: float = 0.0


@dataclass
class Config:
    hours: float = 24.0
    step_s: float = 1.0
    min_elevation_deg: float = 0.0  # visibility mask
    inc_min_deg: float = 15.0  # SAR incidence window (configurable assumption)
    inc_max_deg: float = 50.0


@dataclass
class Geometry:
    t: np.ndarray  # seconds since start
    elevation: np.ndarray
    incidence: np.ndarray
    off_nadir: np.ndarray
    range_m: np.ndarray
    range_rate: np.ndarray  # m/s
    side: np.ndarray  # +1 right looking, -1 left
    sat_lat: np.ndarray
    sat_lon: np.ndarray
    sat_alt: np.ndarray


@dataclass
class Pass:
    sat_name: str
    norad: int
    i0: int
    i1: int  # inclusive
    imaging: list[tuple[int, int]] = field(default_factory=list)


def propagate_ecef(tle: TLE, start: datetime, hours: float, step_s: float):
    """Propagate with SGP4. Returns (t_s, jd, r_ecef[m], v_ecef[m/s])."""
    sat = Satrec.twoline2rv(tle.line1, tle.line2)
    start = start.astimezone(timezone.utc)
    n = int(hours * 3600.0 / step_s) + 1
    t = np.arange(n) * step_s
    jd0, fr0 = jday(start.year, start.month, start.day, start.hour, start.minute,
                    start.second + start.microsecond * 1e-6)
    jd = np.full(n, jd0)
    fr = fr0 + t / 86400.0
    err, r, v = sat.sgp4_array(jd, fr)
    if np.any(err != 0):
        bad = np.flatnonzero(err)[0]
        raise RuntimeError(f"SGP4 error code {err[bad]} for {tle.name} at t={t[bad]:.0f}s")
    r_ecef, v_ecef = geo.teme_to_ecef(r * 1e3, v * 1e3, jd + fr)
    return t, r_ecef, v_ecef


def geometry(t, r_ecef, v_ecef, target: Target) -> Geometry:
    """Look geometry and radial velocity of the satellite as seen from the target."""
    tgt = geo.geodetic_to_ecef(target.lat, target.lon, target.h)
    up_t = geo.up_vector(target.lat, target.lon)

    los = r_ecef - tgt  # target -> satellite
    rng = np.linalg.norm(los, axis=1)
    u = los / rng[:, None]

    cos_inc = np.clip(u @ up_t, -1.0, 1.0)
    incidence = np.degrees(np.arccos(cos_inc))
    elevation = 90.0 - incidence
    rdot = np.einsum("ij,ij->i", u, v_ecef)

    lat_s, lon_s, alt_s = geo.ecef_to_geodetic(r_ecef)
    up_s = geo.up_vector(lat_s, lon_s)
    cos_off = np.clip(np.einsum("ij,ij->i", -u, -up_s), -1.0, 1.0)
    off_nadir = np.degrees(np.arccos(cos_off))

    vhat = v_ecef / np.linalg.norm(v_ecef, axis=1)[:, None]
    right = np.cross(vhat, up_s)  # right-hand side of the ground track
    side = np.where(np.einsum("ij,ij->i", -u, right) >= 0.0, 1, -1)

    return Geometry(t, elevation, incidence, off_nadir, rng, rdot, side, lat_s, lon_s, alt_s)


def doppler_hz(range_rate, freq_hz):
    """Two-way Doppler shift of the carrier for a stationary target."""
    return -2.0 * freq_hz * np.asarray(range_rate) / geo.C


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """Inclusive index ranges where mask is True."""
    if not mask.any():
        return []
    d = np.diff(mask.astype(np.int8), prepend=0, append=0)
    starts = np.flatnonzero(d == 1)
    ends = np.flatnonzero(d == -1) - 1
    return list(zip(starts.tolist(), ends.tolist()))


def find_passes(g: Geometry, tle: TLE, cfg: Config) -> list[Pass]:
    """Contiguous visibility intervals (elevation > mask) with their SAR-eligible sub-windows."""
    passes = []
    for i0, i1 in _runs(g.elevation > cfg.min_elevation_deg):
        seg = slice(i0, i1 + 1)
        ok = (g.incidence[seg] >= cfg.inc_min_deg) & (g.incidence[seg] <= cfg.inc_max_deg)
        imaging = [(i0 + a, i0 + b) for a, b in _runs(ok)]
        passes.append(Pass(tle.name, tle.norad, i0, i1, imaging))
    return passes


def range_accel_at(g: Geometry, idx: int) -> float:
    """d(rdot)/dt [m/s^2] by central difference around sample idx."""
    lo, hi = max(idx - 1, 0), min(idx + 1, len(g.t) - 1)
    return float((g.range_rate[hi] - g.range_rate[lo]) / (g.t[hi] - g.t[lo]))


def utc_iso(start: datetime, seconds: float) -> str:
    return (start + timedelta(seconds=float(seconds))).astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
