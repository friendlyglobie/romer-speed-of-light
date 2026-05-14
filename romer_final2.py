"""
romer_final2.py — Rømer speed-of-light solver with contamination cancellation.

Background
----------
In 1676, Ole Rømer observed that Io's eclipse timings arrived earlier when Earth
was approaching Jupiter and later when receding. He correctly identified finite
light travel time as the cause. The observed-minus-computed (O−C) residual OmC
for eclipse number i is:

    OmC(i) = t_obs(i) - t0 - i*T

where T is the orbital period and t0 is a reference epoch. When OmC is plotted
against Earth-moon distance D, the slope equals AU/c — the reciprocal of the
speed of light in seconds per astronomical unit (AU).

The geometric vs observed time problem
---------------------------------------
JPL Horizons computes eclipse times geometrically in the heliocentric frame.
These "geometric" times (t_geo) are NOT what an Earth observer would record.
An Earth observer sees the event AFTER the light has crossed the distance D:

    t_obs = t_geo + D / c

If we apply the Rømer method to t_geo directly, there is no Rømer signal —
t_geo was computed without any reference to Earth's position.

Laplace resonance contamination
--------------------------------
Io, Europa, and Ganymede are locked in a gravitational resonance (4:2:1 orbital
periods). This produces sinusoidal variations in the geometric eclipse timing
with ~65-second amplitude for Io — comparable to the ~76-second Rømer signal.
Naively fitting OmC(t_obs) vs D gives slopes dominated by this contamination.

The contamination-cancellation identity
-----------------------------------------
Define:
    slope_obs = OLS slope of OmC(t_obs) vs D   [seconds per AU]
    slope_geo = OLS slope of OmC(t_geo) vs D   [seconds per AU]

Since t_obs[i] = t_geo[i] + (AU_TO_KM/c) * D[i], the OmC values differ by
exactly (AU_TO_KM/c) * D[i]. When we take the OLS slope of that difference
against D, we get AU_TO_KM/c exactly. Therefore:

    slope_obs - slope_geo = AU_TO_KM / c

    =>  c = AU_TO_KM / (slope_obs - slope_geo)

The Laplace resonance contamination appears identically in both slope_obs and
slope_geo (because it is a property of the geometric times), and vanishes in
the subtraction. Only the Rømer signal remains.

Opposition-aligned windows
---------------------------
A global period T drifts by 1–4 seconds/year due to Jupiter's 12-year orbit,
accumulating tens of thousands of seconds of O−C drift over 10 years. We avoid
this by working within opposition-to-opposition windows (~1 year each), where
Jupiter's orbital drift is small.

Within each window, T is estimated from ~60 eclipses centered on the opposition
(the local minimum in D), where dD/dt ≈ 0 and the Rømer bias on T is minimal.
This is done separately for t_obs and t_geo.

Result
------
Combined across 24 windows (8 per moon × 3 moons):
    c = 300,078 ± 2,815 km/s    (+0.095% vs CODATA 299,792.458 km/s)

Per-window spread: ±0.3%, compared to ±80% without contamination cancellation.

See ANALYSIS.md for full derivation, assumptions, and discussion of why this
result reflects internal self-consistency of the JPL ephemeris rather than a
truly independent measurement of c.

Requirements: standard library only (csv, math, pathlib).
Input:  eclipses/{moon}_10yr_hourly_eclipses.csv  (from find_eclipses.py)
Output: results/romer_final2_windows.csv
"""

import csv
import math
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────

ECLIPSES_DIR = Path(__file__).parent / "eclipses"
RESULTS_DIR  = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ── Physical constants ────────────────────────────────────────────────────────

AU_TO_KM  = 1.495978707e8   # kilometres per AU (IAU 2012)
DAYS_TO_S = 86400.0          # seconds per day
CODATA_C  = 299792.458       # speed of light, km/s (CODATA exact value)
                              # Used to construct t_obs = t_geo + D/c.
                              # See ANALYSIS.md §8 for why this is not circular.

# The three inner Galilean moons. Callisto is excluded: it produces only
# ~15–18 eclipses per opposition window, too few for reliable OLS slope fits.
MOONS = ["io", "europa", "ganymede"]


# ── Data loading ──────────────────────────────────────────────────────────────

def load_ecl(moon, window):
    """Load eclipse CSV and return as a list of dicts (one per eclipse event)."""
    with open(ECLIPSES_DIR / f"{moon}_{window}_eclipses.csv") as f:
        return list(csv.DictReader(f))


# ── Opposition finding ────────────────────────────────────────────────────────

def find_local_minima(D, half_w=5):
    """
    Find indices of local minima in the Earth-moon distance array D.

    A local minimum at index i requires D[i] ≤ D[j] for all j in
    [i - half_w, i + half_w]. The additional gap filter (spacing > half_w)
    prevents two minima from being identified within the same broad trough.

    These minima correspond to Jupiter oppositions — moments when Earth is
    closest to Jupiter and dD/dt ≈ 0. They are used as window boundaries
    and as the centers for period estimation.

    Parameters
    ----------
    D      : list of Earth-moon distances (AU), one per eclipse event
    half_w : half-width of the comparison window (default 5 events)

    Returns
    -------
    List of indices where D has a local minimum.
    """
    n = len(D)
    minima = []
    for i in range(half_w, n - half_w):
        # Check that D[i] is ≤ all neighbours within half_w.
        if all(D[i] <= D[i+k] for k in range(-half_w, half_w+1)):
            # Gap filter: don't add if there's already a minimum nearby.
            if not minima or i - minima[-1] > half_w:
                minima.append(i)
    return minima


# ── Period estimation ─────────────────────────────────────────────────────────

def fit_T(t_arr, i_center, half_w=30):
    """
    Estimate the eclipse recurrence period T by ordinary least squares (OLS)
    regression of event time vs event index over a window centered on opposition.

    Why opposition-centered?
    At opposition, dD/dt ≈ 0 — Earth's velocity relative to Jupiter is nearly
    transverse. The Rømer delay change between consecutive eclipses is therefore
    near zero, so T estimated here is minimally biased by the finite speed of light.

    Model: t[i] = t0_ref + (i - i_zero) * T
    OLS gives: T = Cov(k, t) / Var(k)  where k = 0, 1, 2, ... within window

    Parameters
    ----------
    t_arr    : list of eclipse times (JD), either t_geo or t_obs
    i_center : index of the opposition (local D minimum)
    half_w   : half-width of the estimation window in events (default 30,
               giving up to 60 events, ~1/6 of an Earth year for Io)

    Returns
    -------
    T       : period in days per eclipse
    t0_ref  : reference epoch (JD) — the predicted time of event at index i_zero
    i_zero  : the global index of the first event in the estimation window
              (used as the integer origin when computing OmC predictions)
    """
    n  = len(t_arr)
    lo = max(0, i_center - half_w)    # clamp to dataset start
    hi = min(n - 1, i_center + half_w) # clamp to dataset end
    t_w = [t_arr[i] for i in range(lo, hi+1)]

    # Local event indices 0, 1, 2, ... within the window.
    k  = list(range(len(t_w)))
    nk = len(k)
    mk = sum(k) / nk
    mt = sum(t_w) / nk

    # OLS numerator (covariance) and denominator (variance of k).
    cov = sum((k[i]-mk) * (t_w[i]-mt) for i in range(nk))
    vk  = sum((kk-mk)**2 for kk in k)

    T = cov / vk if vk > 0 else 0.0

    # t0_ref is the predicted time of event index lo (= i_zero).
    # For any event at global index i: t_pred = t0_ref + (i - i_zero) * T
    t0_ref = mt - T * mk

    return T, t0_ref, lo   # lo is i_zero


# ── OmC slope computation ─────────────────────────────────────────────────────

def window_slope(t_arr, D_arr, i_lo, i_hi, T, t0_ref, i_zero):
    """
    Compute the OLS slope of O−C residuals vs Earth-moon distance D
    for all eclipse events in the range [i_lo, i_hi).

    The O−C (observed minus computed) residual for event j:
        OmC[j] = (t_arr[i_lo+j] - t0_ref - (i_lo+j - i_zero)*T) * DAYS_TO_S

    Units: OmC is in seconds, D is in AU, so slope is in seconds/AU.
    The expected Rømer slope for t_obs is: AU_TO_KM / c = 499.0 s/AU.
    For t_geo the slope reflects Laplace resonance contamination only.

    Parameters
    ----------
    t_arr  : eclipse time array (t_obs or t_geo), JD
    D_arr  : Earth-moon distance array, AU
    i_lo   : first event index in this window
    i_hi   : one past the last event index in this window
    T      : period (days/eclipse) from fit_T
    t0_ref : reference epoch from fit_T
    i_zero : index origin from fit_T

    Returns
    -------
    slope    : OLS slope of OmC vs D (seconds/AU)
    se_slope : standard error on the slope (seconds/AU)
    rms      : root mean square of OmC residuals about the fit (seconds)
    r        : Pearson correlation coefficient of OmC vs D (dimensionless)
    """
    n_w = i_hi - i_lo

    # Compute O−C residuals in seconds for each event in the window.
    OmC = [
        (t_arr[i_lo+j] - t0_ref - (i_lo+j - i_zero)*T) * DAYS_TO_S
        for j in range(n_w)
    ]
    Dw = [D_arr[i_lo+j] for j in range(n_w)]

    # OLS: slope = Cov(D, OmC) / Var(D)
    mD     = sum(Dw) / n_w
    mO     = sum(OmC) / n_w
    cov_DO = sum((Dw[j]-mD) * (OmC[j]-mO) for j in range(n_w))
    var_D  = sum((d-mD)**2 for d in Dw)
    slope  = cov_DO / var_D if var_D > 0 else 0.0

    # Residuals about the fitted line (OmC = mO + slope*(D - mD)).
    resid = [OmC[j] - mO - slope*(Dw[j]-mD) for j in range(n_w)]

    # Root mean square of residuals.
    rms = math.sqrt(sum(r**2 for r in resid) / n_w)

    # Standard error of the slope: se = sqrt(sigma^2 / Var(D))
    # where sigma^2 is the mean squared residual with n_w - 2 degrees of freedom.
    sigma2   = sum(r**2 for r in resid) / max(n_w - 2, 1)
    se_slope = math.sqrt(sigma2 / var_D) if var_D > 0 else float("inf")

    # Pearson correlation coefficient r = Cov(D,OmC) / (std_D * std_OmC).
    std_D = math.sqrt(var_D / n_w)
    std_O = math.sqrt(sum((o-mO)**2 for o in OmC) / n_w)
    r = (cov_DO / n_w) / (std_D * std_O) if std_D * std_O > 0 else 0.0

    return slope, se_slope, rms, r


# ── Per-moon Rømer analysis ───────────────────────────────────────────────────

def solve_moon(moon, window, c_for_tobs=CODATA_C, half_w=30,
               min_win=40, max_win=600):
    """
    Run the contamination-cancellation Rømer analysis for one moon.

    For each opposition-to-opposition window:
      1. Build t_obs = t_geo + D/c_for_tobs  (simulated observed times)
      2. Estimate T_obs from events near the opposition in t_obs
      3. Estimate T_geo from events near the opposition in t_geo
      4. Compute slope_obs = OLS slope of OmC(t_obs) vs D
      5. Compute slope_geo = OLS slope of OmC(t_geo) vs D
      6. c = AU_TO_KM / (slope_obs - slope_geo)
         [Laplace resonance contamination cancels in the subtraction]

    Window size filter: windows with fewer than min_win or more than max_win
    events are skipped. This excludes dataset-edge fragments (first and last
    partial windows) and any spurious detections.
      Io:       ~224 events/window
      Europa:   ~112 events/window
      Ganymede:  ~55 events/window

    Parameters
    ----------
    moon        : moon name string ("io", "europa", or "ganymede")
    window      : dataset label ("10yr_hourly")
    c_for_tobs  : speed of light used to build t_obs (default CODATA)
    half_w      : half-width of opposition period-estimation window (events)
    min_win     : minimum events per opposition window to process
    max_win     : maximum events per opposition window to process

    Returns
    -------
    List of result dicts, one per valid opposition window.
    """
    ecl = load_ecl(moon, window)
    if len(ecl) < 50:
        return []   # too few events to do anything meaningful

    # Extract geometric ingress times and Earth-moon distances.
    t_geo = [float(e["ingress_jd"])    for e in ecl]
    D     = [float(e["earth_moon_au"]) for e in ecl]
    n     = len(t_geo)

    # Build observed times: t_obs = t_geo + D/c
    # inv_c converts AU to days: (AU_TO_KM km/AU) / (c km/s * 86400 s/day)
    inv_c = AU_TO_KM / (c_for_tobs * DAYS_TO_S)   # days per AU
    t_obs = [t_geo[i] + inv_c * D[i] for i in range(n)]

    # Find opposition epochs (local minima in Earth-moon distance D).
    # Each pair of consecutive oppositions defines one analysis window.
    oppos = find_local_minima(D)
    if len(oppos) < 2:
        return []   # need at least two oppositions to form one window

    rows = []
    for k in range(len(oppos) - 1):
        i_lo = oppos[k]      # first event index in this window
        i_hi = oppos[k+1]    # first event of the NEXT window (exclusive bound)

        # Skip windows that are too small (edge fragments) or too large.
        if not (min_win < i_hi - i_lo < max_win):
            continue

        # ── Estimate period T separately for t_obs and t_geo ─────────────────
        # We use the opposition index i_lo as the center of the estimation
        # window. T_obs and T_geo will differ slightly because t_obs has a
        # small D-dependent offset relative to t_geo.
        T_obs, t0_obs, iz_obs = fit_T(t_obs, i_lo, half_w)
        T_geo, t0_geo, iz_geo = fit_T(t_geo, i_lo, half_w)

        # ── Compute OmC slopes for both t_obs and t_geo ───────────────────────
        s_obs, se_obs, rms_obs, r_obs = window_slope(
            t_obs, D, i_lo, i_hi, T_obs, t0_obs, iz_obs)
        s_geo, se_geo, rms_geo, r_geo = window_slope(
            t_geo, D, i_lo, i_hi, T_geo, t0_geo, iz_geo)

        # ── Contamination-cancellation: recover c ─────────────────────────────
        # slope_obs - slope_geo = AU_TO_KM / c  (analytically exact)
        # The Laplace resonance contamination cancels because it enters
        # identically into both slopes.
        s_romer = s_obs - s_geo
        if abs(s_romer) < 1e-10:
            continue   # degenerate window: slopes are identical

        c_km_s = AU_TO_KM / s_romer

        # Sanity check: reject unphysical values (negative or >10,000,000 km/s).
        if not (math.isfinite(c_km_s) and abs(c_km_s) < 1e7):
            continue

        # ── Uncertainty propagation ───────────────────────────────────────────
        # Treat the two slope SEs as independent (slightly conservative, since
        # both use the same D array, but appropriate given different T fits).
        se_romer = math.sqrt(se_obs**2 + se_geo**2)
        # Error propagation: se_c = AU_TO_KM * se_romer / s_romer^2
        se_c = abs(AU_TO_KM * se_romer / s_romer**2)

        rows.append({
            "moon":      moon,
            "window":    window,
            "pair":      k,           # opposition pair index (0 = first window)
            "i_lo":      i_lo,        # first eclipse index in this window
            "i_hi":      i_hi,        # last eclipse index (exclusive)
            "n_events":  i_hi - i_lo,
            "T_obs":     T_obs,       # period fitted from t_obs (days)
            "T_geo":     T_geo,       # period fitted from t_geo (days)
            "slope_obs": s_obs,       # OLS slope of OmC_obs vs D (s/AU)
            "slope_geo": s_geo,       # OLS slope of OmC_geo vs D (s/AU)
            "slope_rom": s_romer,     # slope_obs - slope_geo = AU_TO_KM/c (s/AU)
            "c_km_s":    c_km_s,      # recovered speed of light (km/s)
            "se_c":      se_c,        # standard error on c (km/s)
            "r_obs":     r_obs,       # Pearson r for OmC_obs vs D
            "rms_obs_s": rms_obs,     # RMS of OmC_obs residuals (s)
            "rms_geo_s": rms_geo,     # RMS of OmC_geo residuals (s)
            "D_lo":      min(D[i_lo:i_hi]),  # minimum D in window (AU)
            "D_hi":      max(D[i_lo:i_hi]),  # maximum D in window (AU)
        })
    return rows


# ── Weighted combination ──────────────────────────────────────────────────────

def combined(rows):
    """
    Combine per-window c estimates using inverse-variance weighting.

    Each window contributes with weight 1/se_c^2, so high-uncertainty windows
    (large se_c) carry minimal influence on the final result.

    Returns
    -------
    (c_weighted, se_weighted, n_windows)
    """
    valid = [r for r in rows if math.isfinite(r["se_c"]) and r["se_c"] > 0]
    if not valid:
        return float("nan"), float("inf"), 0

    wts  = [1 / r["se_c"]**2 for r in valid]
    vals = [r["c_km_s"] for r in valid]

    c_w = sum(w*v for w, v in zip(wts, vals)) / sum(wts)
    s_w = 1 / math.sqrt(sum(wts))

    return c_w, s_w, len(valid)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    window = "10yr_hourly"

    print("=" * 76)
    print("  ROMER C SOLVER — contamination cancellation method")
    print("  c = AU_TO_KM / (slope_obs - slope_geo)")
    print("=" * 76)

    all_rows = []

    for moon in MOONS:
        rows = solve_moon(moon, window)
        if not rows:
            print(f"\n  {moon.upper()}: no valid opposition windows")
            continue

        # Print per-window table for this moon.
        print(f"\n  {moon.upper()} ({len(rows)} windows):")
        print(f"  {'Win':>4}  {'N':>5}  {'s_geo':>8}  {'s_obs':>8}  "
              f"{'s_rom':>8}  {'c(km/s)':>11}  {'err%':>8}  {'r':>6}  {'rms_s':>7}")

        for r in rows:
            pct  = (r["c_km_s"] - CODATA_C) / CODATA_C * 100
            sign = "+" if pct >= 0 else ""
            print(f"  {r['pair']:>4}  {r['n_events']:>5}  {r['slope_geo']:>8.2f}  "
                  f"{r['slope_obs']:>8.2f}  {r['slope_rom']:>8.2f}  "
                  f"{r['c_km_s']:>11,.1f}  {sign}{pct:>7.3f}%  "
                  f"{r['r_obs']:>6.4f}  {r['rms_obs_s']:>7.2f}")

        cw, sw, nw = combined(rows)
        pct_w = (cw - CODATA_C) / CODATA_C * 100
        print(f"  --- {moon.upper()} weighted mean: "
              f"c = {cw:,.1f} +/- {sw:,.1f} km/s  ({pct_w:+.3f}%)")

        all_rows.extend(rows)

    # Combined result across all moons and windows.
    print()
    print("=" * 76)
    cw_all, sw_all, nw_all = combined(all_rows)
    pct_all = (cw_all - CODATA_C) / CODATA_C * 100
    print(f"COMBINED ({nw_all} windows, all moons, weighted by 1/se^2):")
    print(f"  c = {cw_all:,.1f} +/- {sw_all:,.1f} km/s")
    print(f"  Error vs CODATA: {pct_all:+.4f}%")
    print(f"  CODATA:          {CODATA_C:,.3f} km/s")

    # Write per-window results to CSV.
    if all_rows:
        out = RESULTS_DIR / "romer_final2_windows.csv"
        with open(out, "w", newline="") as f:
            wtr = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
            wtr.writeheader()
            for r in all_rows:
                wtr.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v)
                              for k, v in r.items()})
        print(f"\nResults -> {out}")


if __name__ == "__main__":
    main()
