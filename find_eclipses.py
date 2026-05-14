"""
find_eclipses.py — Detect Galilean moon umbral ingress times from JPL vectors.

Overview
--------
For each Galilean moon we scan hourly heliocentric position vectors and find
every moment the moon enters Jupiter's umbral shadow cone. These geometric
ingress times are the raw material for the Rømer speed-of-light analysis.

Shadow geometry (all coordinates Jupiter-centered)
---------------------------------------------------
Working in Jupiter-centered coordinates simplifies the shadow cone calculation.

    r_moon = pos_moon - pos_jupiter     (moon relative to Jupiter, AU)
    r_sun  = -pos_jupiter               (Sun relative to Jupiter, AU;
                                         the Sun sits at the heliocentric origin)

Jupiter's umbra is a cone extending away from the Sun. At axial distance d
from Jupiter along the shadow axis, the umbra radius is:

    R_umbra(d) = R_jup - d * (R_sun - R_jup) / D_sun_jup

This narrows linearly because the Sun is larger than Jupiter (R_sun > R_jup).
The cone eventually closes to a point (the umbra tip) at:

    d_tip = R_jup * D_sun_jup / (R_sun - R_jup)  ≈ 7.3 AU from Jupiter

All three inner Galilean moons orbit well within this distance.

A moon is inside the umbra when both:
  1. It is on the anti-solar side of Jupiter: dot(r_moon, shadow_axis) > 0
  2. Its perpendicular distance from the shadow axis < R_umbra at its depth

Shadow depth scalar:
    depth = R_umbra(d) - r_perp
    depth > 0  =>  inside umbra (eclipsed)
    depth < 0  =>  outside

Ingress detection
-----------------
We evaluate depth at every hourly sample. A sign change from ≤0 to >0
brackets an ingress within a 1-hour interval. Brent's method applied to a
Hermite cubic spline of the position vectors then narrows the ingress time
to < 10 ms precision.

Why Hermite splines?
JPL Horizons provides both position and velocity at each sample. A Hermite
cubic spline uses both, matching position AND velocity at the endpoints.
This gives physically realistic interpolation between samples and is much
more accurate than linear interpolation for smooth orbital trajectories.

Important: ingress_jd is a GEOMETRIC time
------------------------------------------
The ingress_jd stored in the output CSV is the moment the moon crosses the
umbra boundary in the heliocentric inertial frame. Earth plays no role in this
calculation. An observer on Earth would see the event LATER, after light has
traveled the Earth-moon distance D at speed c:

    t_obs = t_geo + D / c    (D in AU, c in AU/s)

The Rømer solver (romer_final2.py) constructs t_obs from t_geo and D.

Output columns per eclipse event
---------------------------------
  moon              — moon name
  ingress_jd        — geometric ingress time (JD TDB), <10 ms precision
  ingress_cal       — human-readable calendar equivalent
  egress_jd         — geometric egress time (JD TDB)
  egress_cal        — human-readable calendar equivalent
  duration_min      — eclipse duration in minutes
  mid_jd            — midpoint JD
  earth_moon_au     — Earth-to-moon distance at ingress (AU)
  light_travel_min  — light travel time Earth-moon at ingress (minutes)
  ingress_prec_s    — Brent convergence precision in seconds (~0.009 s)
  egress_prec_s     — same for egress

Physical constants
------------------
  R_sun = 0.00465047 AU  (696,000 km solar radius)
  R_jup = 0.000477895 AU (71,492 km Jupiter equatorial radius)

Requirements: standard library only (csv, math, pathlib).
Input:  vectors/{moon}_{window}.csv, vectors/jupiter_{window}.csv,
        vectors/earth_{window}.csv  (produced by fetch_vectors.py)
Output: eclipses/{moon}_{window}_eclipses.csv
"""

import csv
import math
from pathlib import Path

# ── Physical constants ────────────────────────────────────────────────────────

R_SUN_AU  = 0.00465047     # solar radius in AU (696,000 km)
R_JUP_AU  = 0.000477895    # Jupiter equatorial radius in AU (71,492 km)
AU_TO_KM  = 1.495978707e8  # kilometres per AU (IAU 2012)
C_KM_S    = 299792.458     # speed of light in km/s (CODATA exact)
C_AU_MIN  = C_KM_S * 60 / AU_TO_KM   # speed of light in AU/minute

# Brent's method stops when the bracket is narrower than this.
# 1e-7 days = 8.64 ms — well below any physical uncertainty in the geometry.
BRENT_TOL_DAYS = 1e-7

# ── Paths ─────────────────────────────────────────────────────────────────────

BASE = Path(__file__).parent / "vectors"    # input: state vector CSVs
OUT  = Path(__file__).parent / "eclipses"   # output: eclipse event CSVs
OUT.mkdir(exist_ok=True)

# Moons and time windows to process.
# The Rømer analysis uses only the three inner moons and "10yr_hourly".
# Callisto has too few eclipses per opposition window for reliable slope fits.
MOONS   = ["io", "europa", "ganymede", "callisto"]
WINDOWS = ["1yr_hourly", "10yr_hourly"]


# ── Vector arithmetic (plain tuples, no numpy dependency) ─────────────────────

def sub(a, b):
    """Vector subtraction: a - b."""
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def dot(a, b):
    """Vector dot product."""
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def norm(a):
    """Euclidean magnitude of a 3-vector."""
    return math.sqrt(dot(a, a))

def vadd(a, b):
    """Vector addition."""
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])

def vscale(a, s):
    """Scalar multiplication of a 3-vector."""
    return (a[0]*s, a[1]*s, a[2]*s)


# ── Shadow geometry ───────────────────────────────────────────────────────────

def shadow_depth(r_moon_jup, r_sun_jup) -> float:
    """
    Compute the shadow depth scalar for a moon in Jupiter's umbra.

    Returns a value > 0 if the moon is inside the umbra, < 0 if outside.

    Parameters
    ----------
    r_moon_jup : (x, y, z) tuple — moon position relative to Jupiter, AU
    r_sun_jup  : (x, y, z) tuple — Sun position relative to Jupiter, AU
                 (= -pos_jupiter in heliocentric coords, since Sun is at origin)
    """
    # Distance from Jupiter to Sun; used to scale the umbra cone.
    D = norm(r_sun_jup)

    # Unit vector pointing from Jupiter toward the Sun.
    sun_hat = vscale(r_sun_jup, 1.0/D)

    # Shadow axis points in the OPPOSITE direction: away from the Sun.
    axis = vscale(sun_hat, -1.0)

    # Project the moon's position onto the shadow axis.
    # d > 0 means the moon is on the anti-solar (shadow) side of Jupiter.
    d = dot(r_moon_jup, axis)
    if d <= 0:
        # Moon is on the sunlit side — definitely not in shadow.
        return -1.0

    # Perpendicular distance of the moon from the shadow axis.
    proj   = vscale(axis, d)                   # component along axis
    r_perp = norm(sub(r_moon_jup, proj))        # component perpendicular to axis

    # Umbra cone radius at axial depth d.
    # The formula is linear: R_jup at d=0, narrowing toward zero at the tip.
    R_umbra = R_JUP_AU - d * (R_SUN_AU - R_JUP_AU) / D
    if R_umbra <= 0:
        # Moon is beyond the tip of the umbra cone.
        return -1.0

    # Positive = inside umbra, negative = outside.
    return R_umbra - r_perp


# ── Hermite cubic spline interpolation ───────────────────────────────────────

def cubic_interp_vec(t, t0, t1, p0, p1, v0, v1):
    """
    Interpolate a 3-vector position using a Hermite cubic spline.

    Uses both endpoint positions (p0, p1) and endpoint velocities (v0, v1)
    to construct a cubic that matches position AND velocity at t0 and t1.
    This is the natural choice when velocities are available from JPL Horizons
    and gives physically realistic trajectories between hourly samples.

    Parameters
    ----------
    t          : query time (JD)
    t0, t1     : bracket endpoint times (JD)
    p0, p1     : 3-vector positions at t0, t1 (AU)
    v0, v1     : 3-vector velocities at t0, t1 (AU/day)

    Returns
    -------
    Interpolated 3-vector position at time t (AU).
    """
    h  = t1 - t0            # interval width in days
    s  = (t - t0) / h       # normalised parameter in [0, 1]
    s2 = s * s
    s3 = s2 * s

    # Hermite basis functions. These are the four cubic polynomials that
    # satisfy: h00(0)=1, h00(1)=0; h10(0)=0 (deriv=1 at 0, 0 at 1); etc.
    h00 =  2*s3 - 3*s2 + 1   # blends toward p0 as s → 0
    h10 =    s3 - 2*s2 + s   # handles velocity at t0 (scaled by h below)
    h01 = -2*s3 + 3*s2       # blends toward p1 as s → 1
    h11 =    s3 -   s2       # handles velocity at t1 (scaled by h below)

    # The h* factor on velocity terms converts from AU/day to the right scale
    # for the normalised parameter s.
    return (
        h00*p0[0] + h10*h*v0[0] + h01*p1[0] + h11*h*v1[0],
        h00*p0[1] + h10*h*v0[1] + h01*p1[1] + h11*h*v1[1],
        h00*p0[2] + h10*h*v0[2] + h01*p1[2] + h11*h*v1[2],
    )


# ── Brent's method root finder ────────────────────────────────────────────────

def brent(f, xa, xb, fa, fb, tol=BRENT_TOL_DAYS, maxiter=60):
    """
    Find the root of f(t) = 0 in the bracket [xa, xb].

    Brent's method combines three strategies in order of preference:
      1. Inverse quadratic interpolation — cubic-rate convergence near root
      2. Secant method — superlinear convergence
      3. Bisection fallback — guaranteed bracket halving

    It always maintains a valid bracket [xa, xb] where f changes sign,
    so convergence is guaranteed regardless of function shape.

    Parameters
    ----------
    f       : callable, the function to find the root of
    xa, xb  : bracket endpoints (must bracket a root: fa * fb < 0)
    fa, fb  : f(xa) and f(xb), pre-computed to save one function call
    tol     : stop when bracket width < tol (default 1e-7 days = 8.6 ms)
    maxiter : maximum iterations before returning best estimate

    Returns
    -------
    (root, iterations) — root estimate and number of iterations used.
    """
    # Convention: xb is always the current best estimate (smallest |f|).
    if abs(fa) < abs(fb):
        xa, xb = xb, xa
        fa, fb = fb, fa

    xc, fc = xa, fa     # third point, starts at xa
    mflag = True        # True = last step was bisection or IQI
    s = xb
    d = 0.0

    for i in range(maxiter):
        if abs(xb - xa) < tol:
            break

        if fa != fc and fb != fc:
            # Inverse quadratic interpolation through three distinct points.
            s = (xa*fb*fc / ((fa-fb)*(fa-fc))
               + xb*fa*fc / ((fb-fa)*(fb-fc))
               + xc*fa*fb / ((fc-fa)*(fc-fb)))
        else:
            # Secant step through the two best points.
            s = xb - fb * (xb - xa) / (fb - fa)

        # Fall back to bisection if the interpolated step is out of bounds
        # or not shrinking the bracket fast enough.
        cond1 = not ((3*xa + xb)/4 < s < xb or xb < s < (3*xa + xb)/4)
        cond2 = mflag     and abs(s - xb) >= abs(xb - xc) / 2
        cond3 = not mflag and abs(s - xb) >= abs(xc - d) / 2
        cond4 = mflag     and abs(xb - xc) < tol
        cond5 = not mflag and abs(xc - d)  < tol

        if cond1 or cond2 or cond3 or cond4 or cond5:
            s = (xa + xb) / 2   # bisect
            mflag = True
        else:
            mflag = False

        fs = f(s)
        d, xc, fc = xc, xb, fb     # shift history

        # Update bracket: keep whichever side contains the sign change.
        if fa * fs < 0:
            xb, fb = s, fs
        else:
            xa, fa = s, fs

        # Keep xb as the better estimate.
        if abs(fa) < abs(fb):
            xa, xb = xb, xa
            fa, fb = fb, fa

    return xb, i + 1


# ── Data loading ──────────────────────────────────────────────────────────────

def load_body(path: Path):
    """
    Load a JPL Horizons vector CSV into parallel Python lists.

    Returns
    -------
    jds : list of float  — Julian Dates (TDB)
    pos : list of (x,y,z) tuples — heliocentric position (AU)
    vel : list of (vx,vy,vz) tuples — heliocentric velocity (AU/day)
    """
    jds, pos, vel = [], [], []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            jds.append(float(row["jd_tdb"]))
            pos.append((float(row["x_au"]), float(row["y_au"]), float(row["z_au"])))
            vel.append((float(row["vx_au_d"]), float(row["vy_au_d"]), float(row["vz_au_d"])))
    return jds, pos, vel


def jd_to_cal_approx(jd: float) -> str:
    """
    Convert a Julian Date to an approximate Gregorian calendar string.

    Uses the standard astronomical algorithm (Meeus, "Astronomical Algorithms").
    Accurate to about 1 second for dates in the modern era. The output timezone
    is TDB (Barycentric Dynamical Time), which differs from UTC by at most ~2 ms.

    Example output: "A.D. 2022-03-15 14:32:07.41 TDB"
    """
    # JD 2451545.0 = J2000.0 = 2000-Jan-01 12:00 TT
    z = int(jd + 0.5)
    if z < 2299161:
        a = z   # Julian calendar
    else:
        # Gregorian calendar correction
        alpha = int((z - 1867216.25) / 36524.25)
        a = z + 1 + alpha - alpha // 4
    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)

    day_frac = (jd + 0.5) - int(jd + 0.5)
    day   = b - d - int(30.6001 * e)
    month = e - 1 if e < 14 else e - 13
    year  = c - 4716 if month > 2 else c - 4715

    total_s = day_frac * 86400
    hh = int(total_s // 3600)
    mm = int((total_s % 3600) // 60)
    ss = total_s % 60

    return f"A.D. {year}-{month:02d}-{day:02d} {hh:02d}:{mm:02d}:{ss:05.2f} TDB"


# ── Main eclipse detection ────────────────────────────────────────────────────

def process(moon: str, window: str):
    """
    Detect all umbral ingress and egress events for one moon in one time window.

    Algorithm
    ---------
    1. Load position+velocity arrays for the moon, Jupiter, and Earth.
    2. Compute shadow_depth at every hourly sample.
    3. Scan for sign changes: depth ≤0 → >0 is an ingress; >0 → ≤0 is an egress.
    4. For each bracketed sign change, build a Hermite spline and apply Brent's
       method to find the exact crossing time to <10 ms precision.
    5. Record the Earth-moon distance (by spline-interpolating Earth's position
       to the exact ingress time) for use in the Rømer analysis.
    6. Write results to eclipses/{moon}_{window}_eclipses.csv.
    """
    # Load state vectors. All three bodies must share the same time grid.
    moon_jds, moon_pos, moon_vel = load_body(BASE / f"{moon}_{window}.csv")
    jup_jds,  jup_pos,  jup_vel  = load_body(BASE / f"jupiter_{window}.csv")
    erth_jds, erth_pos, erth_vel = load_body(BASE / f"earth_{window}.csv")

    n = len(moon_jds)
    assert len(jup_jds) == n == len(erth_jds), "Row count mismatch between bodies"

    # ── Step 1: pre-compute shadow depth at every hourly sample ──────────────
    # This is cheap (pure arithmetic) and lets us scan for sign changes in O(n).
    depths = []
    for i in range(n):
        # Convert heliocentric positions to Jupiter-centered coordinates.
        r_mj = sub(moon_pos[i], jup_pos[i])    # moon relative to Jupiter
        r_sj = vscale(jup_pos[i], -1.0)        # Sun relative to Jupiter
                                                # (Sun is at heliocentric origin,
                                                #  so Sun_jup = 0 - jup = -jup)
        depths.append(shadow_depth(r_mj, r_sj))

    # ── Step 2: factory for the refined depth function ───────────────────────
    # Given a bracketing pair of sample indices, returns a closure f(t) that
    # spline-interpolates both moon and Jupiter positions to any time t within
    # the bracket and evaluates shadow_depth at that interpolated state.
    # Brent's method calls this closure to narrow the crossing time.
    def make_depth_fn(ia, ib):
        t0, t1 = moon_jds[ia], moon_jds[ib]
        mp0, mp1 = moon_pos[ia], moon_pos[ib]
        mv0, mv1 = moon_vel[ia], moon_vel[ib]
        jp0, jp1 = jup_pos[ia],  jup_pos[ib]
        jv0, jv1 = jup_vel[ia],  jup_vel[ib]

        def f(t):
            # Spline-interpolate moon and Jupiter positions to time t.
            pm = cubic_interp_vec(t, t0, t1, mp0, mp1, mv0, mv1)
            pj = cubic_interp_vec(t, t0, t1, jp0, jp1, jv0, jv1)
            # Re-compute Jupiter-centered coordinates from interpolated positions.
            r_mj = sub(pm, pj)
            r_sj = vscale(pj, -1.0)
            return shadow_depth(r_mj, r_sj)

        return f

    def earth_moon_dist_at(jd: float, ia: int, ib: int) -> float:
        """
        Compute Earth-to-moon distance (AU) at a precise JD by spline interpolation.

        Used to record the Earth-moon range D at the exact geometric ingress time.
        D is later used by the Rømer solver: t_obs = t_geo + D/c.
        """
        t0, t1 = moon_jds[ia], moon_jds[ib]
        pm = cubic_interp_vec(jd, t0, t1, moon_pos[ia], moon_pos[ib],
                              moon_vel[ia], moon_vel[ib])
        pe = cubic_interp_vec(jd, t0, t1, erth_pos[ia], erth_pos[ib],
                              erth_vel[ia], erth_vel[ib])
        return norm(sub(pm, pe))

    # ── Step 3: scan for sign changes and refine crossing times ──────────────
    events = []
    in_eclipse  = depths[0] > 0    # are we inside the umbra at t=0?
    ingress_jd  = None              # refined ingress time, set when bracket found
    ingress_ia  = None              # index of the sample before ingress bracket

    for i in range(1, n):
        pd, cd = depths[i-1], depths[i]    # previous and current depth

        if not in_eclipse and cd > 0 and pd <= 0:
            # ── Ingress: depth crossed from ≤0 to >0 ─────────────────────────
            # Build the spline depth function for this 1-hour bracket and
            # apply Brent's method to find the exact crossing time.
            f = make_depth_fn(i-1, i)
            ingress_jd, ingress_iters = brent(f, moon_jds[i-1], moon_jds[i], pd, cd)
            ingress_ia = i - 1
            in_eclipse = True

        elif in_eclipse and cd <= 0 and pd > 0:
            # ── Egress: depth crossed from >0 to ≤0 ──────────────────────────
            f = make_depth_fn(i-1, i)
            egress_jd, egress_iters = brent(f, moon_jds[i-1], moon_jds[i], pd, cd)
            in_eclipse = False

            # Only record the event if we have a valid matching ingress.
            if ingress_jd is not None:
                mid_jd       = (ingress_jd + egress_jd) / 2
                duration_min = (egress_jd - ingress_jd) * 24 * 60

                # Earth-moon distance at the precise ingress moment.
                # This D value is what the Rømer solver uses to build t_obs.
                dist_au = earth_moon_dist_at(ingress_jd, ingress_ia, ingress_ia + 1)
                lt_min  = dist_au / C_AU_MIN   # light travel time in minutes

                # Precision estimate: the Brent tolerance converted to seconds.
                ingress_prec_s = BRENT_TOL_DAYS * 86400   # ~0.009 s
                egress_prec_s  = BRENT_TOL_DAYS * 86400

                events.append({
                    "moon":             moon,
                    "ingress_jd":       f"{ingress_jd:.9f}",
                    "ingress_cal":      jd_to_cal_approx(ingress_jd),
                    "egress_jd":        f"{egress_jd:.9f}",
                    "egress_cal":       jd_to_cal_approx(egress_jd),
                    "duration_min":     f"{duration_min:.4f}",
                    "mid_jd":           f"{mid_jd:.9f}",
                    "earth_moon_au":    f"{dist_au:.8f}",   # key input for Rømer
                    "light_travel_min": f"{lt_min:.6f}",
                    "ingress_prec_s":   f"{ingress_prec_s:.4f}",
                    "egress_prec_s":    f"{egress_prec_s:.4f}",
                })
            # Reset ingress state for the next eclipse.
            ingress_jd = None
            ingress_ia = None

    # ── Step 4: write output CSV ──────────────────────────────────────────────
    out_path = OUT / f"{moon}_{window}_eclipses.csv"
    if events:
        with open(out_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(events[0].keys()))
            w.writeheader()
            w.writerows(events)
        print(f"  {moon:10s} {window:14s}: {len(events):4d} eclipses -> {out_path.name}")
    else:
        print(f"  {moon:10s} {window:14s}: no eclipses found")


def main():
    """Process all moons for all time windows."""
    for window in WINDOWS:
        print(f"\n=== {window} ===")
        for moon in MOONS:
            process(moon, window)


if __name__ == "__main__":
    main()
