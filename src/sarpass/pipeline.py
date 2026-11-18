"""End-to-end: TLEs -> passes over a target -> JSON consumed by the web page."""

from __future__ import annotations

from datetime import datetime, timezone

import numpy as np

from . import sar
from .sar import Config, Mission, Target
from .tle import TLE

EXPORT_STEP_S = 2  # sampling of the curves stored in the JSON


def _r(a, nd):
    return np.round(np.asarray(a, dtype=float), nd).tolist()


def build_output(
    tles: list[TLE],
    mission: Mission,
    target: Target,
    cfg: Config,
    start: datetime,
    source: str,
    max_tle_age_days: float = 30.0,
) -> dict:
    start = start.astimezone(timezone.utc)
    stride = max(int(round(EXPORT_STEP_S / cfg.step_s)), 1)
    passes_out: list[dict] = []
    sats_out: list[dict] = []
    skipped: list[str] = []
    visible_not_imageable = 0

    for tle in tles:
        age = (start - tle.epoch).total_seconds() / 86400.0
        if abs(age) > max_tle_age_days:
            skipped.append(f"{tle.name} (TLE epoch {age:+.0f} d)")
            continue
        try:
            t, r, v = sar.propagate_ecef(tle, start, cfg.hours, cfg.step_s)
        except RuntimeError as exc:  # decayed / bad elements
            skipped.append(f"{tle.name} ({exc})")
            continue
        g = sar.geometry(t, r, v, target)
        passes = sar.find_passes(g, tle, cfg)
        sats_out.append(
            {"name": tle.name, "norad": tle.norad, "tle_epoch": tle.epoch.strftime("%Y-%m-%dT%H:%M:%SZ"),
             "tle_age_days": round(age, 2), "passes_visible": len(passes),
             "passes_imageable": sum(1 for p in passes if p.imaging)}
        )

        for p in passes:
            if not p.imaging:
                visible_not_imageable += 1
                continue
            seg = slice(p.i0, p.i1 + 1)
            i_tca = p.i0 + int(np.argmin(g.range_m[seg]))
            sub = np.arange(p.i0, p.i1 + 1, stride)
            t0 = g.t[p.i0]
            passes_out.append(
                {
                    "sat": tle.name,
                    "norad": tle.norad,
                    "start_utc": sar.utc_iso(start, g.t[p.i0]),
                    "end_utc": sar.utc_iso(start, g.t[p.i1]),
                    "tca_utc": sar.utc_iso(start, g.t[i_tca]),
                    "duration_s": int(g.t[p.i1] - g.t[p.i0]),
                    "truncated": bool(p.i0 == 0 or p.i1 == len(g.t) - 1),
                    "t_tca_s": float(g.t[i_tca] - t0),
                    "max_elevation_deg": round(float(g.elevation[seg].max()), 2),
                    "min_range_km": round(float(g.range_m[i_tca] / 1e3), 1),
                    "incidence_tca_deg": round(float(g.incidence[i_tca]), 2),
                    "off_nadir_tca_deg": round(float(g.off_nadir[i_tca]), 2),
                    "sat_alt_km": round(float(g.sat_alt[i_tca] / 1e3), 1),
                    "side": "right" if g.side[i_tca] > 0 else "left",
                    "range_accel_ms2": round(sar.range_accel_at(g, i_tca), 4),
                    "imaging_windows_s": [[float(g.t[a] - t0), float(g.t[b] - t0)] for a, b in p.imaging],
                    # curves
                    "t_s": _r(g.t[sub] - t0, 0),
                    "range_km": _r(g.range_m[sub] / 1e3, 1),
                    "range_rate_ms": _r(g.range_rate[sub], 1),
                    "incidence_deg": _r(g.incidence[sub], 2),
                    "sat_lat": _r(g.sat_lat[sub], 3),
                    "sat_lon": _r(g.sat_lon[sub], 3),
                }
            )

    passes_out.sort(key=lambda d: d["tca_utc"])
    for i, p in enumerate(passes_out):
        p["id"] = i

    return {
        "schema": 1,
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "window_start_utc": start.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "window_hours": cfg.hours,
        "source": source,  # "celestrak" | "demo" | "file"
        "mission": {
            "key": mission.key,
            "label": mission.label,
            "band": mission.band,
            "freq_hz": mission.freq_hz,
            "wavelength_cm": round(mission.wavelength_m * 100.0, 3),
        },
        "target": {"name": target.name, "lat": target.lat, "lon": target.lon},
        "params": {
            "inc_min_deg": cfg.inc_min_deg,
            "inc_max_deg": cfg.inc_max_deg,
            "min_elevation_deg": cfg.min_elevation_deg,
            "step_s": cfg.step_s,
        },
        "satellites": sats_out,
        "skipped": skipped,
        "passes_visible_not_imageable": visible_not_imageable,
        "passes": passes_out,
    }
