"""TLE utilities: parsing, checksum, and a generator for synthetic (demo) TLEs."""

from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime, timezone

MU = 398600.4418  # km^3/s^2
R_E_KM = 6378.137


@dataclass
class TLE:
    name: str
    line1: str
    line2: str

    @property
    def norad(self) -> int:
        return int(self.line1[2:7])

    @property
    def epoch(self) -> datetime:
        yy = int(self.line1[18:20])
        year = 2000 + yy if yy < 57 else 1900 + yy
        doy = float(self.line1[20:32])
        return datetime(year, 1, 1, tzinfo=timezone.utc) + _days(doy - 1.0)


def _days(d: float):
    from datetime import timedelta

    return timedelta(days=d)


def checksum(line: str) -> int:
    """Modulo-10 TLE checksum of the first 68 characters."""
    total = 0
    for ch in line[:68]:
        if ch.isdigit():
            total += int(ch)
        elif ch == "-":
            total += 1
    return total % 10


def parse_tle_text(text: str) -> list[TLE]:
    """Parse 3-line (name + 2 lines) or 2-line TLE text. Bad-checksum sets are skipped."""
    lines = [ln.rstrip() for ln in text.splitlines() if ln.strip()]
    out: list[TLE] = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("1 ") and i + 1 < len(lines) and lines[i + 1].startswith("2 "):
            name, l1, l2 = f"NORAD {lines[i][2:7].strip()}", lines[i], lines[i + 1]
            i += 2
        elif (
            i + 2 < len(lines)
            and lines[i + 1].startswith("1 ")
            and lines[i + 2].startswith("2 ")
        ):
            name = lines[i].strip()
            if name.startswith("0 "):  # CelesTrak 3LE prefix
                name = name[2:].strip()
            l1, l2 = lines[i + 1], lines[i + 2]
            i += 3
        else:
            i += 1
            continue
        if len(l1) >= 69 and len(l2) >= 69 and checksum(l1) == int(l1[68]) and checksum(l2) == int(l2[68]):
            out.append(TLE(name, l1[:69], l2[:69]))
    return out


def make_tle(
    name: str,
    satnum: int,
    epoch: datetime,
    altitude_km: float,
    inclination_deg: float,
    raan_deg: float,
    mean_anomaly_deg: float,
    ecc: float = 0.0001,
    argp_deg: float = 90.0,
) -> TLE:
    """Build a syntactically valid TLE from circular-orbit parameters (SYNTHETIC data)."""
    a = R_E_KM + altitude_km
    n_rad_s = math.sqrt(MU / a**3)
    rev_per_day = n_rad_s * 86400.0 / (2.0 * math.pi)
    epoch = epoch.astimezone(timezone.utc)
    doy = (epoch - datetime(epoch.year, 1, 1, tzinfo=timezone.utc)).total_seconds() / 86400.0 + 1.0
    l1 = (
        f"1 {satnum:05d}U 00000A   {epoch.year % 100:02d}{doy:012.8f}"
        f"  .00000000  00000-0  00000-0 0  999"
    )
    l1 += str(checksum(l1 + "0"))
    l2 = (
        f"2 {satnum:05d} {inclination_deg:8.4f} {raan_deg % 360:8.4f} "
        f"{int(round(ecc * 1e7)):07d} {argp_deg % 360:8.4f} {mean_anomaly_deg % 360:8.4f} "
        f"{rev_per_day:11.8f}{1:5d}"
    )
    l2 += str(checksum(l2 + "0"))
    assert len(l1) == 69 and len(l2) == 69, (len(l1), len(l2))
    return TLE(name, l1, l2)
