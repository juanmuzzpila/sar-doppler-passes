"""Command line interface:  python -m sarpass --mission iceye --out docs/data/passes.json"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from . import celestrak
from .demo import demo_constellation
from .pipeline import build_output
from .sar import MISSIONS, Config, Target, doppler_hz
from .tle import parse_tle_text


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="sarpass", description=__doc__)
    ap.add_argument("--mission", choices=sorted(MISSIONS), default="iceye")
    ap.add_argument("--lat", type=float, default=39.4699, help="target latitude [deg] (default: Valencia)")
    ap.add_argument("--lon", type=float, default=-0.3763, help="target longitude [deg] (default: Valencia)")
    ap.add_argument("--name", default="Valencia, ES")
    ap.add_argument("--hours", type=float, default=24.0)
    ap.add_argument("--step", type=float, default=1.0, help="propagation step [s]")
    ap.add_argument("--inc-min", type=float, default=15.0)
    ap.add_argument("--inc-max", type=float, default=50.0)
    ap.add_argument("--min-elev", type=float, default=0.0)
    ap.add_argument("--start", help="window start, ISO-8601 UTC (default: now)")
    ap.add_argument("--demo", action="store_true", help="use a synthetic constellation (no network)")
    ap.add_argument("--tle-file", type=Path, help="use TLEs from this file instead of CelesTrak")
    ap.add_argument("--save-tles", type=Path, help="after a successful CelesTrak download, also save the TLEs to this file")
    ap.add_argument("--cache-dir", type=Path, default=Path(".cache"))
    ap.add_argument("--max-tle-age-days", type=float, default=30.0)
    ap.add_argument("--out", type=Path, default=Path("docs/data/passes.json"))
    args = ap.parse_args(argv)

    mission = MISSIONS[args.mission]
    start = (
        datetime.fromisoformat(args.start.replace("Z", "+00:00")).astimezone(timezone.utc)
        if args.start
        else datetime.now(timezone.utc).replace(microsecond=0)
    )

    if args.demo:
        tles, source = demo_constellation(start), "demo"
    elif args.tle_file:
        tles, source = parse_tle_text(args.tle_file.read_text()), "file"
    else:
        try:
            tles, origin = celestrak.fetch_by_name(mission.celestrak_query, args.cache_dir)
        except Exception as exc:  # noqa: BLE001
            print(f"error: could not get TLEs from CelesTrak ({exc}).\n"
                  f"       Re-run with --demo or --tle-file <file>.", file=sys.stderr)
            return 2
        source = "celestrak"
        if args.save_tles:
            args.save_tles.parent.mkdir(parents=True, exist_ok=True)
            args.save_tles.write_text("".join(f"{t.name}\n{t.line1}\n{t.line2}\n" for t in tles))
            print(f"saved {len(tles)} TLEs -> {args.save_tles}")
        print(f"TLEs: {len(tles)} for NAME={mission.celestrak_query} ({origin})")

    if not tles:
        print("error: no TLEs to process", file=sys.stderr)
        return 2

    cfg = Config(hours=args.hours, step_s=args.step, min_elevation_deg=args.min_elev,
                 inc_min_deg=args.inc_min, inc_max_deg=args.inc_max)
    target = Target(args.name, args.lat, args.lon)
    out = build_output(tles, mission, target, cfg, start, source, args.max_tle_age_days)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, separators=(",", ":")))

    print(f"Target: {target.name} ({target.lat:.4f}, {target.lon:.4f}) | "
          f"{mission.label} {mission.band}-band {mission.freq_hz/1e9:.3f} GHz | source={source}")
    print(f"Window: {out['window_start_utc']} + {cfg.hours:g} h | imaging incidence {cfg.inc_min_deg:g}-{cfg.inc_max_deg:g} deg")
    print(f"{len(out['passes'])} SAR-eligible passes ({out['passes_visible_not_imageable']} more visible but outside incidence range)")
    for p in out["passes"]:
        rr = np.array(p["range_rate_ms"])
        fd = doppler_hz(rr, mission.freq_hz) / 1e3
        print(f"  {p['tca_utc']}  {p['sat'][:28]:<28}  {p['side']:<5}  inc@TCA {p['incidence_tca_deg']:5.1f}  "
              f"max el {p['max_elevation_deg']:5.1f}  Doppler {fd.min():+7.1f}..{fd.max():+7.1f} kHz")
    if out["skipped"]:
        print(f"skipped {len(out['skipped'])} object(s): " + "; ".join(out["skipped"][:5]))
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
