"""Synthetic demo constellation (NOT real satellites).

Used when CelesTrak is unreachable (offline, CI sandbox, first look at the page).
Six sun-synchronous satellites at ~570 km with staggered planes/phasing so the
demo shows a few passes per day. Everything produced from these is labelled
`source: "demo"` in the output JSON and flagged in the web UI.
"""

from __future__ import annotations

import math
from datetime import datetime

from .tle import MU, R_E_KM, TLE, make_tle

J2 = 1.08263e-3
SSO_NODAL_RATE = 2.0 * math.pi / (365.2422 * 86400.0)  # rad/s


def sso_inclination_deg(altitude_km: float) -> float:
    """Inclination for which J2 nodal precession matches the Sun (circular orbit)."""
    a = R_E_KM + altitude_km
    n = math.sqrt(MU / a**3)
    cos_i = -2.0 * SSO_NODAL_RATE * a**2 / (3.0 * J2 * n * R_E_KM**2)
    return math.degrees(math.acos(cos_i))


def demo_constellation(epoch: datetime, n_sats: int = 6, altitude_km: float = 570.0) -> list[TLE]:
    inc = sso_inclination_deg(altitude_km)
    tles = []
    for k in range(n_sats):
        tles.append(
            make_tle(
                name=f"DEMO-SAR-{k + 1} (synthetic)",
                satnum=99001 + k,
                epoch=epoch,
                altitude_km=altitude_km,
                inclination_deg=inc,
                raan_deg=40.0 + k * 55.0,
                mean_anomaly_deg=k * 137.0,
            )
        )
    return tles
