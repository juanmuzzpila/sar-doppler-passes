"""CelesTrak downloader with an on-disk cache (be polite to their service).

CelesTrak asks clients not to re-download the same data more often than every
~2 hours, so results are cached and re-used within that window.
"""

from __future__ import annotations

import time
import urllib.parse
import urllib.request
from pathlib import Path

from .tle import TLE, parse_tle_text

URL = "https://celestrak.org/NORAD/elements/gp.php"
CACHE_TTL_S = 2 * 3600
UA = "sar-doppler-passes (educational project; github.com)"


def fetch_by_name(query: str, cache_dir: Path = Path(".cache")) -> tuple[list[TLE], str]:
    """Return (TLEs, origin) where origin is 'network', 'cache' or 'stale-cache'."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    safe = "".join(c if c.isalnum() else "_" for c in query)
    path = cache_dir / f"celestrak_{safe}.tle"

    if path.exists() and time.time() - path.stat().st_mtime < CACHE_TTL_S:
        return parse_tle_text(path.read_text()), "cache"

    url = f"{URL}?{urllib.parse.urlencode({'NAME': query, 'FORMAT': 'tle'})}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=30) as resp:
            text = resp.read().decode("utf-8", "replace")
        tles = parse_tle_text(text)
        if not tles:
            raise RuntimeError(f"CelesTrak returned no usable TLEs for NAME={query!r}: {text[:80]!r}")
        path.write_text(text)
        return tles, "network"
    except Exception:
        if path.exists():  # fall back to whatever we have
            return parse_tle_text(path.read_text()), "stale-cache"
        raise
