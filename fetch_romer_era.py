"""
fetch_romer_era.py
==================
Download heliocentric state vectors from JPL Horizons for the three bodies
needed to compute Io eclipse times during the Rømer era (1671–1677).

WHY HELIOCENTRIC VECTORS?
  Eclipse detection requires knowing where Io sits inside Jupiter's shadow cone.
  The shadow cone is defined by the Sun–Jupiter–Io geometry, so we need positions
  of all three bodies in a common inertial frame (the Solar System barycentre /
  heliocentric ecliptic J2000).  We also need Earth's position to compute the
  light-travel time from Io to the observer.

BODIES FETCHED:
  Earth   (Horizons ID 399)  — observer position; needed for light-travel time
  Jupiter (Horizons ID 5)    — shadow-cone apex; Io position is relative to this
  Io      (Horizons ID 501)  — the moon whose eclipses Rømer timed

TIME WINDOW:
  1671-01-01 to 1677-01-01, step = 1 hour.
  Io's orbital period is ~1.769 days, so 1-hour sampling gives ~42 points per
  orbit — enough to bracket every ingress/egress within a single hour, which
  is then refined to millisecond precision by Hermite interpolation + Brent root-finding
  in find_eclipses_romer_era.py.

VECTOR TABLE TYPE 2 (position + velocity):
  We request both position (x, y, z) and velocity (vx, vy, vz) in AU and AU/day.
  The velocity is essential for the Hermite cubic spline interpolation used in
  find_eclipses_romer_era.py — with both position and velocity at each endpoint,
  a cubic Hermite polynomial uniquely matches position AND derivative at both ends,
  giving far better accuracy than position-only linear interpolation.

OUTPUT:
  vectors/earth_romer_era.csv   — 52,609 hourly rows, 1671–1677
  vectors/jupiter_romer_era.csv
  vectors/io_romer_era.csv

Run this script once; it skips bodies whose output file already exists.
Then run find_eclipses_romer_era.py to detect eclipse times.
"""

import time
import csv
import sys
import warnings
import requests
from pathlib import Path

# JPL returns self-signed HTTPS; suppress the warning so output is clean.
warnings.filterwarnings("ignore", message="Unverified HTTPS")

HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
OUT_DIR      = Path(__file__).parent / "vectors"
OUT_DIR.mkdir(exist_ok=True)

# Horizons body IDs for the three bodies we need.
BODIES = {
    "earth":   "399",
    "jupiter": "5",
    "io":      "501",
}

# Six years covering Rømer's observing campaign at the Paris Observatoire Royal.
WINDOW = {
    "label": "romer_era",
    "start": "1671-01-01",
    "stop":  "1677-01-01",
    "step":  "1 h",
}

# Horizons API parameters common to all three requests.
#
# CENTER "@10"   = Sun (body 10), so vectors are heliocentric.
# REF_PLANE      = ECLIPTIC — ecliptic plane of J2000; consistent with how
#                  JPL's planetary ephemerides are defined.
# VEC_TABLE 2    = position + velocity (needed for Hermite interpolation).
# VEC_CORR NONE  = geometric vectors, no aberration correction.  We handle
#                  light-travel time explicitly in the analysis script.
# OUT_UNITS AU-D = AU for distance, AU/day for velocity.
# CSV_FORMAT YES = machine-readable CSV output between $$SOE/$$EOE markers.
BASE_PARAMS = {
    "format":      "json",
    "EPHEM_TYPE":  "VECTORS",
    "MAKE_EPHEM":  "YES",
    "OBJ_DATA":    "NO",
    "CENTER":      "@10",
    "REF_PLANE":   "ECLIPTIC",
    "VEC_TABLE":   "2",
    "VEC_CORR":    "NONE",
    "OUT_UNITS":   "AU-D",
    "CSV_FORMAT":  "YES",
    "CAL_FORMAT":  "CAL",
}


def fetch(body_id, start, stop, step):
    """
    Query JPL Horizons for one body and return a list of row dicts.

    Each row has keys: jd_tdb, cal_date, x_au, y_au, z_au, vx_au_d, vy_au_d, vz_au_d.
    The Horizons response wraps the data table between $$SOE and $$EOE markers;
    we extract that block and parse it as CSV.
    """
    params = BASE_PARAMS | {
        "COMMAND":    f"'{body_id}'",
        "START_TIME": f"'{start}'",
        "STOP_TIME":  f"'{stop}'",
        "STEP_SIZE":  f"'{step}'",
    }
    resp = requests.get(HORIZONS_URL, params=params, timeout=120, verify=False)
    data = resp.json()

    if "error" in data:
        raise RuntimeError(f"Horizons error: {data['error']}")

    raw = data.get("result", "")
    soe = raw.find("$$SOE")
    eoe = raw.find("$$EOE")
    if soe == -1 or eoe == -1:
        raise RuntimeError("Could not find $$SOE/$$EOE markers in response")

    block = raw[soe + 5 : eoe].strip()
    rows  = []
    for line in block.splitlines():
        line  = line.strip()
        if not line:
            continue
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 8:
            continue
        try:
            rows.append({
                "jd_tdb":   parts[0],   # Julian Date, TDB timescale
                "cal_date": parts[1],   # calendar date string (for human reading)
                "x_au":     parts[2],   # heliocentric x position (AU)
                "y_au":     parts[3],   # heliocentric y position (AU)
                "z_au":     parts[4],   # heliocentric z position (AU)
                "vx_au_d":  parts[5],   # x velocity (AU/day)
                "vy_au_d":  parts[6],   # y velocity (AU/day)
                "vz_au_d":  parts[7],   # z velocity (AU/day)
            })
        except IndexError:
            continue
    return rows


def save_csv(rows, path):
    if not rows:
        print(f"  WARNING: no rows for {path.name}", flush=True)
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"  Saved {len(rows):,} rows -> {path.name}", flush=True)


def main():
    win = WINDOW
    print(f"Fetching Romer era vectors ({win['start']} to {win['stop']}, step={win['step']})\n",
          flush=True)
    for name, body_id in BODIES.items():
        out_path = OUT_DIR / f"{name}_{win['label']}.csv"
        if out_path.exists():
            print(f"  {name}: already exists, skipping", flush=True)
            continue
        print(f"  Fetching {name} ({body_id})...", end=" ", flush=True)
        try:
            rows = fetch(body_id, win["start"], win["stop"], win["step"])
            print(f"{len(rows):,} rows", flush=True)
            save_csv(rows, out_path)
        except Exception as e:
            print(f"FAILED: {e}", file=sys.stderr, flush=True)
        time.sleep(1)   # be polite to the Horizons API
    print("\nDone. Run find_eclipses_romer_era.py to detect eclipse times.", flush=True)


if __name__ == "__main__":
    main()
