# Rømer's Method: Speed of Light from Galilean Moon Eclipses

A from-scratch Python replication of Ole Rømer's 1676 measurement of the speed
of light, using modern ephemeris data from the Jet Propulsion Laboratory (JPL)
Horizons system instead of a telescope.

## Result

```
c = 300,078 ± 2,815 km/s   (+0.095% vs CODATA 299,792.458 km/s)
```

Recovered across 24 independent opposition windows spanning three Galilean
moons (Io, Europa, Ganymede) over a 10-year baseline (2016–2026).

## How it works

Rømer noticed in 1676 that Io's eclipses arrive early when Earth is approaching
Jupiter and late when receding. The delay is the light travel time across the
changing Earth–Jupiter distance — roughly ±8 minutes over one Earth year.

The central challenge is that the Laplace orbital resonance among Io, Europa,
and Ganymede produces timing variations of comparable amplitude (~65 seconds)
that would otherwise swamp the Rømer signal. This analysis cancels that
contamination exactly using the identity:

```
c  =  AU / (slope_obs − slope_geo)
```

where `slope_obs` is the ordinary least squares (OLS) slope of the observed-minus-computed
(O−C) residual vs Earth–moon distance using Earth-observed eclipse times, and
`slope_geo` is the same slope using the geometric (JPL heliocentric frame)
eclipse times. The Laplace contamination appears identically in both slopes and
vanishes in the difference, leaving only the Rømer signal.

For a full derivation, all assumptions, and a discussion of why this analysis
demonstrates internal self-consistency rather than a truly independent
measurement of c, see [ANALYSIS.md](ANALYSIS.md).

## Reproduce

**Requirements**: Python 3.9+, standard library only (no third-party packages).

```bash
# Step 1 — fetch JPL Horizons state vectors (requires internet, ~5 minutes)
python fetch_vectors.py

# Step 2 — detect eclipse ingress times from shadow geometry
python find_eclipses.py

# Step 3 — run the Rømer solver
python romer_final2.py
```

Step 1 writes hourly heliocentric position and velocity vectors for Earth,
Jupiter, Io, Europa, Ganymede, and Callisto into `vectors/` (gitignored —
large, but fully reproducible). Steps 2 and 3 use only files already in this
repository if you skip Step 1 (the `eclipses/` CSVs are pre-computed and
included).

## Files

| File | Purpose |
|------|---------|
| `fetch_vectors.py` | Queries JPL Horizons API for 10-year hourly state vectors |
| `find_eclipses.py` | Detects umbral ingress times via Hermite spline + Brent's method |
| `romer_final2.py` | Contamination-cancellation Rømer solver |
| `eclipses/` | Pre-computed eclipse timing CSVs (one per moon) |
| `results/romer_final2_windows.csv` | Per-window c estimates, slopes, and uncertainties |
| `ANALYSIS.md` | Full derivation, method, assumptions, and results |

## Per-moon summary

| Moon     | Windows | c (km/s) | Error vs CODATA |
|----------|---------|----------|-----------------|
| Io       | 8       | 299,921  | +0.043%         |
| Europa   | 8       | 299,949  | +0.052%         |
| Ganymede | 8       | 300,203  | +0.137%         |
| Combined | 24      | 300,078  | +0.095%         |

## A note on independence

The JPL ephemeris positions used here were derived from decades of radar
ranging, spacecraft transponder timing, and Very Long Baseline Interferometry
(VLBI) measurements — all of which assumed the CODATA value of c in their data
reduction. The Earth–moon distances in `eclipses/` therefore already embed c
implicitly. What this analysis demonstrates is that the contamination-cancellation
method correctly recovers the c that went in, which is a non-trivial
verification of the method. A genuinely independent measurement would require
raw telescope timing records reduced without any prior assumption about c —
exactly what Rømer had in 1676.
