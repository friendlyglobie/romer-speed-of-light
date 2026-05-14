"""
fetch_vectors.py — Download heliocentric state vectors from JPL Horizons.

Fetches position and velocity vectors for Earth, Jupiter, and the four Galilean
moons (Io, Europa, Ganymede, Callisto) from NASA's JPL Horizons web API.

All vectors are in the heliocentric ecliptic J2000 reference frame: positions
in astronomical units (AU), velocities in AU/day, times in Julian Date (JD)
Barycentric Dynamical Time (TDB).

The primary dataset used by the Rømer solver is "10yr_hourly": ten years at
1-hour resolution. The other windows are provided for reference.

Output: CSV files written to ./vectors/, one file per body per time window.
Each CSV has columns: jd_tdb, cal_date, x_au, y_au, z_au, vx_au_d, vy_au_d, vz_au_d

The vectors/ directory is excluded from the repository (see .gitignore) because
the files are large (~400 MB total for 10yr_hourly). Run this script once to
regenerate them locally before running find_eclipses.py.

Requirements: pip install requests
"""

import time
import csv
import sys
import warnings
import requests
from pathlib import Path

# JPL Horizons sometimes presents a self-signed certificate; suppress the warning.
warnings.filterwarnings("ignore", message="Unverified HTTPS")

# ── Constants and configuration ──────────────────────────────────────────────

HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"

# Output directory for raw vector CSV files (gitignored — regenerate with this script).
OUT_DIR = Path(__file__).parent / "vectors"
OUT_DIR.mkdir(exist_ok=True)

# JPL Horizons body identifiers.
# Earth = 399, Jupiter barycenter = 5 (includes all satellites).
# Galilean moons: Io = 501, Europa = 502, Ganymede = 503, Callisto = 504.
BODIES = {
    "earth":    "399",
    "jupiter":  "5",
    "io":       "501",
    "europa":   "502",
    "ganymede": "503",
    "callisto": "504",
}

# Time windows to fetch. The Rømer analysis uses "10yr_hourly" exclusively.
# "10yr" (daily) is a lightweight reference dataset.
# "1yr_hourly" covers a single opposition and is too short for the solver
# (which requires at least two oppositions to define a window).
WINDOWS = [
    {
        "label": "10yr",
        "start": "2016-01-01",
        "stop":  "2026-01-01",
        "step":  "1 d",           # daily — small file, useful for inspection
    },
    {
        "label": "1yr_hourly",
        "start": "2024-09-01",
        "stop":  "2025-09-01",
        "step":  "1 h",           # hourly over one year
    },
    {
        "label": "10yr_hourly",
        "start": "2016-01-01",
        "stop":  "2026-01-01",
        "step":  "1 h",           # hourly over ten years — primary dataset
    },
]

# Parameters sent to every Horizons API request.
# CENTER=@10 means heliocentric (Sun = body 10).
# REF_PLANE=ECLIPTIC gives ecliptic J2000 coordinates.
# VEC_TABLE=2 returns X Y Z VX VY VZ (positions + velocities).
# OUT_UNITS=AU-D gives positions in AU and velocities in AU/day.
BASE_PARAMS = {
    "format":      "json",
    "EPHEM_TYPE":  "VECTORS",
    "MAKE_EPHEM":  "YES",
    "OBJ_DATA":    "NO",
    "CENTER":      "@10",          # heliocentric origin
    "REF_PLANE":   "ECLIPTIC",
    "VEC_TABLE":   "2",            # position + velocity
    "VEC_CORR":    "NONE",         # geometric vectors, no aberration correction
    "OUT_UNITS":   "AU-D",         # AU and AU/day
    "CSV_FORMAT":  "YES",
    "CAL_FORMAT":  "CAL",
}


# ── API fetch ─────────────────────────────────────────────────────────────────

def fetch(body_id: str, start: str, stop: str, step: str) -> list[dict]:
    """
    Query JPL Horizons for one body over one time window.

    Horizons returns a text block between $$SOE (Start Of Ephemeris) and $$EOE
    (End Of Ephemeris) markers. With CSV_FORMAT=YES each line is:
        JDTDB, Calendar_Date(TDB), X, Y, Z, VX, VY, VZ, [extra cols...]

    Returns a list of dicts, one per time step, with keys:
        jd_tdb, cal_date, x_au, y_au, z_au, vx_au_d, vy_au_d, vz_au_d
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

    # Locate the ephemeris data block between the $$SOE and $$EOE sentinel lines.
    soe = raw.find("$$SOE")
    eoe = raw.find("$$EOE")
    if soe == -1 or eoe == -1:
        raise RuntimeError("Could not find $$SOE/$$EOE markers in response")

    block = raw[soe + 5 : eoe].strip()
    rows = []
    for line in block.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = [p.strip() for p in line.split(",")]
        # Need at least 8 fields: JD, date, x, y, z, vx, vy, vz.
        if len(parts) < 8:
            continue
        try:
            rows.append({
                "jd_tdb":  parts[0],   # Julian Date, TDB timescale
                "cal_date": parts[1],  # human-readable calendar date
                "x_au":    parts[2],   # heliocentric X position (AU)
                "y_au":    parts[3],   # heliocentric Y position (AU)
                "z_au":    parts[4],   # heliocentric Z position (AU)
                "vx_au_d": parts[5],   # X velocity (AU/day)
                "vy_au_d": parts[6],   # Y velocity (AU/day)
                "vz_au_d": parts[7],   # Z velocity (AU/day)
            })
        except IndexError:
            continue
    return rows


# ── Output ────────────────────────────────────────────────────────────────────

def save_csv(rows: list[dict], path: Path) -> None:
    """Write a list of row dicts to a CSV file."""
    if not rows:
        print(f"  WARNING: no rows to write to {path.name}")
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"  Saved {len(rows):,} rows -> {path.name}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    """
    Fetch all bodies for all time windows, skipping files that already exist.
    A 1-second pause between requests avoids hammering the Horizons API.
    """
    for win in WINDOWS:
        print(f"\n=== Window: {win['label']} "
              f"({win['start']} to {win['stop']}, step={win['step']}) ===")
        for name, body_id in BODIES.items():
            out_path = OUT_DIR / f"{name}_{win['label']}.csv"
            if out_path.exists():
                print(f"  {name}: already exists, skipping")
                continue
            print(f"  Fetching {name} ({body_id})...", end=" ", flush=True)
            try:
                rows = fetch(body_id, win["start"], win["stop"], win["step"])
                print(f"{len(rows):,} rows")
                save_csv(rows, out_path)
            except Exception as e:
                print(f"FAILED: {e}", file=sys.stderr)
            time.sleep(1)  # be polite to the Horizons API


if __name__ == "__main__":
    main()
