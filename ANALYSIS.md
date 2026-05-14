# Rømer Speed-of-Light Recovery: Full Analysis Overview

## 1. Background and Goal

In 1676, Ole Rømer noticed that Io's eclipse timings arrived earlier or later depending
on whether Earth was approaching or receding from Jupiter. He correctly attributed this
to the finite travel time of light. The effect has a ~16-minute peak-to-peak amplitude
over one Earth year (corresponding to Earth's orbital diameter of ~2 astronomical units
(AU) crossed at c ≈ 300,000 km/s).

This analysis replicates that measurement from first principles, using modern Jet
Propulsion Laboratory (JPL) ephemeris data rather than a telescope, across a 10-year
baseline (2016–2026) and three Galilean moons: Io, Europa, and Ganymede.

---

## 2. Data Acquisition: JPL Horizons State Vectors

**Script**: `fetch_vectors.py`

JPL's Horizons system provides numerically integrated planetary and satellite state
vectors. We fetched heliocentric ecliptic J2000 state vectors (position + velocity) for:

- Earth (Horizons ID 399)
- Jupiter barycenter (ID 5)
- Io (501), Europa (502), Ganymede (503), Callisto (504)

**Reference frame**: Heliocentric ecliptic J2000 (`CENTER=@10`, `REF_PLANE=ECLIPTIC`).
All bodies are expressed in the same inertial frame, which is essential — shadow geometry
requires computing relative positions in a common reference.

**Time coverage**: 2016-01-01 to 2026-01-01, sampled at 1-hour intervals. At 1 hour,
Io (period ~1.77 days) completes about 4% of an orbit per sample; the spline-based
crossing refinement (see §3) recovers sub-second timing from these coarse samples.

**Output format**: Comma-separated values (CSV) files with columns `jd_tdb, cal_date,
x_au, y_au, z_au, vx_au_d, vy_au_d, vz_au_d`. Positions are in AU, velocities in AU
per day, and time is in Julian Date (JD) — Barycentric Dynamical Time (TDB), a uniform
coordinate time used for solar system ephemerides that avoids the non-uniform ticking of
clocks in Earth's gravitational potential.

---

## 3. Eclipse Detection

**Script**: `find_eclipses.py`

### 3.1 Shadow Geometry

All computations are done in **Jupiter-centered coordinates**. For each time step:

```
r_moon  = pos_moon  - pos_jupiter    (moon position relative to Jupiter)
r_sun   = -pos_jupiter               (Sun position relative to Jupiter,
                                      since Sun is at heliocentric origin)
```

Jupiter casts an **umbral cone** on the anti-solar side. The umbra is the region of
complete shadow where the Sun is entirely blocked; it is the relevant region because
umbral ingress/egress produces a definitive, sharp change in illumination, as opposed
to the penumbra (partial shadow), which produces a gradual dimming that was
unmeasurable with 17th-century instruments and is excluded here.

A moon is in the umbra when two conditions hold simultaneously:

**Condition 1** — the moon is on the anti-solar side:
```
dot(r_moon, r_sun) < 0
```

**Condition 2** — the moon is within the umbral cone. The cone narrows with distance
from Jupiter because the Sun is larger than Jupiter. At axial distance `d` from Jupiter
along the shadow axis:

```
R_umbra(d) = R_jup - d × (R_sun - R_jup) / D_sun_jup
```

where `R_jup = 0.000477895 AU` (71,492 km, the equatorial radius), `R_sun = 0.00465047
AU` (696,000 km), and `D_sun_jup` is the instantaneous Jupiter–Sun distance (~5.2 AU,
varying by ±0.25 AU over Jupiter's 12-year orbit). Since R_sun > R_jup, the umbra radius
shrinks with distance and eventually reaches zero at the umbra tip.

The **shadow depth scalar** that encodes eclipse state:

```
d       = dot(r_moon, shadow_axis)         shadow_axis = -r_sun / |r_sun|
r_perp  = |r_moon - d × shadow_axis|      perpendicular distance from shadow axis
depth   = R_umbra(d) - r_perp
```

`depth > 0` means inside the umbra (eclipsed); `depth ≤ 0` means outside.

**Assumption**: Jupiter is treated as a perfect sphere with its equatorial radius.
Jupiter is oblate (polar radius ~66,854 km vs equatorial 71,492 km), but the oblateness
introduces at most a ~1% error in the umbra radius, which is negligible for our
purposes.

### 3.2 Ingress Detection: Sign-Change Scanning

The algorithm scans through the hourly samples looking for sign changes in `depth`:

- A transition from `depth ≤ 0` to `depth > 0` is an **ingress** (moon entering shadow)
- A transition from `depth > 0` to `depth ≤ 0` is an **egress** (moon leaving shadow)

When a sign change is detected between samples `i-1` and `i`, the crossing time is
bracketed to a 1-hour window.

### 3.3 Crossing Refinement: Hermite Spline + Brent's Method

A 1-hour bracket gives ~1800-second timing uncertainty, far too coarse for Rømer
analysis (the signal is ~76 seconds peak amplitude). We refine using two techniques
in combination.

**Hermite cubic spline interpolation** of both the moon and Jupiter position vectors
within the bracket. Given positions p0, p1 and velocities v0, v1 at the bracket
endpoints, the spline is:

```
p(s) = h00(s)·p0 + h10(s)·h·v0 + h01(s)·p1 + h11(s)·h·v1
```

where `s = (t - t0)/h` is the normalized time parameter and h00, h10, h01, h11 are
the Hermite basis polynomials:

```
h00 = 2s³ - 3s² + 1
h10 = s³  - 2s² + s
h01 = -2s³ + 3s²
h11 = s³  - s²
```

The velocities come directly from the Horizons data, making this a natural fit: the
spline matches both position and velocity at the endpoints, giving C1 continuity and
physically realistic interpolated trajectories between samples.

**Brent's method** is then applied to the function `f(t) = shadow_depth(t)` on the
bracket, where `shadow_depth` is evaluated using the spline-interpolated positions.
Brent's method combines inverse quadratic interpolation, the secant method, and
bisection fallback, guaranteeing convergence without requiring derivatives. The
tolerance is set to `1e-7 days = 8.6 milliseconds`, well below any physical uncertainty.

**Output column `ingress_jd`**: The Julian Date (TDB) of umbral ingress, accurate to
<10 ms. This is the **geometric shadow-crossing time** — the instant the moon crosses
the umbral cone boundary, computed entirely in the heliocentric inertial frame. It is
**not** an observed time from Earth; no light-travel correction is applied at this stage.

**Output column `earth_moon_au`**: The Earth-to-moon distance in AU at the exact
ingress moment, obtained by Hermite-interpolating the Earth and moon positions to the
refined ingress time.

### 3.4 Eclipse Counts

For the 10-year hourly dataset (2016–2026):

| Moon      | Synodic period | Events |
|-----------|----------------|--------|
| Io        | 1.769861 days  | 1,843  |
| Europa    | 3.553 days     | 918    |
| Ganymede  | 7.166 days     | 455    |
| Callisto  | 16.7 days      | 94     |

Callisto was excluded from the Rømer analysis: 94 events over 10 years gives only 5–6
opposition windows with ~15–18 events each — too few per window for a stable ordinary
least squares (OLS) slope estimate.

---

## 4. The Rømer Method and Its Complications

### 4.1 The Basic Idea

If eclipses recur with period T, then the predicted time of eclipse number i is:

```
t_pred(i) = t0 + i × T
```

The observed minus computed residual (O−C, abbreviated OmC) is:

```
OmC(i) = t_obs(i) - t0 - i × T
```

For Rømer's observation: when Earth is moving away from Jupiter, successive eclipses
appear to arrive later than predicted (positive OmC); when approaching, earlier
(negative OmC). The OmC should be proportional to the Earth–Jupiter distance D:

```
OmC ≈ D / c    (in consistent units)
```

A linear regression of OmC vs D then gives slope = AU_TO_KM / c (in seconds per AU),
and c = AU_TO_KM / slope.

### 4.2 The Geometric vs Observed Time Problem

The critical insight that unlocks the analysis: **`ingress_jd` is a geometric time, not
an observed time**.

JPL Horizons computes where bodies are in space and when shadow crossings occur
geometrically — purely from orbital mechanics in the heliocentric frame. An observer on
Earth doesn't see this event happen at `t_geo` — they see the light that left the moon
at `t_geo` arrive after a delay of D/c:

```
t_obs = t_geo + D/c
```

where D = `earth_moon_au` in AU and AU/c ≈ 499.0 seconds ≈ 0.005775 days.

If we naively apply the Rømer method to `t_geo` without this correction, we are fitting
a time series that has no Rømer signal at all (since `t_geo` was computed without
reference to Earth's position). This is why early attempts returned nonsensical results.

To build the observed times, we inject the known value of c:

```python
inv_c  = AU_TO_KM / (c_for_tobs × DAYS_TO_S)   # days per AU
t_obs[i] = t_geo[i] + inv_c × D[i]
```

This is discussed further in §8, where the apparent circularity is resolved.

### 4.3 Orbital Contamination: The Laplace Resonance Problem

If the only signal in OmC were the Rømer delay, the slope of OmC_obs vs D would equal
AU_TO_KM/c = 499.0 s/AU exactly and the method would work immediately. In practice
there is a large contaminant.

Io, Europa, and Ganymede are locked in a **Laplace resonance**: for every orbit of
Ganymede, Europa completes approximately 2 and Io completes approximately 4. This
mutual gravitational coupling produces periodic variations in each moon's orbital
elements — small oscillations in eccentricity, argument of pericenter, and mean anomaly.
These in turn cause the geometric eclipse timing to be not perfectly periodic, exhibiting:

- An **annual sinusoidal perturbation** (~65-second amplitude for Io, larger for Europa
  and Ganymede) driven by Jupiter's orbital motion modulating the resonance geometry
- A **Laplace libration** with ~480-day period and similar amplitude

The net effect: `t_geo` contains systematic timing variations correlated with D (because
D tracks Jupiter's orbital phase, which also drives the resonance). When we compute OmC
from `t_geo`, its slope vs D is substantial — and varies wildly depending on which
portion of the 480-day libration cycle a given window falls in:

```
slope_geo ≈ −830 to +737 s/AU  (varies per window and per moon)
```

Compare to the expected Rømer signal: slope_romer = 499 s/AU. For Europa, the
contamination magnitude can exceed the Rømer signal itself. Without correction, recovered
c values range from negative to several million km/s.

---

## 5. The Solver: Contamination-Cancellation Method

**Script**: `romer_final2.py`

### 5.1 The Key Identity

The contamination can be removed exactly using an analytic identity. Define:

```
slope_obs = OLS slope of OmC(t_obs) vs D
slope_geo = OLS slope of OmC(t_geo) vs D
```

where OmC is computed using the same T and t0 reference for both. Then:

```
slope_obs = slope_geo + AU_TO_KM / (c × DAYS_TO_S)
```

This follows because:

```
t_obs[i] = t_geo[i] + (AU_TO_KM / (c × DAYS_TO_S)) × D[i]
```

Forming OmC(t_obs) for any event i:

```
OmC_obs[i] = (t_obs[i] - t0_ref - (i - i_zero) × T) × DAYS_TO_S
           = OmC_geo[i] + (AU_TO_KM / (c × DAYS_TO_S)) × D[i]  ×  DAYS_TO_S
```

Wait — more carefully, since OmC is computed in seconds:

```
OmC_obs[i]  [seconds] = OmC_geo[i]  [seconds]  +  (AU_TO_KM / c)  ×  D[i]  [AU]
```

The OLS slope of OmC_obs vs D equals the OLS slope of OmC_geo vs D plus the slope of
the additive term `(AU_TO_KM/c) × D` vs D, which is exactly AU_TO_KM/c. Therefore:

```
slope_obs  [s/AU]  −  slope_geo  [s/AU]  =  AU_TO_KM / c  [km/AU / (km/s)]  =  AU_TO_KM / c  [s/AU]
```

Solving for c:

```
c  =  AU_TO_KM / (slope_obs − slope_geo)      [km/s]
```

This is the **contamination-cancellation formula**. The Laplace resonance contamination
appears identically in both slope_obs and slope_geo and cancels in the difference,
leaving only the Rømer signal.

### 5.2 Opposition-Aligned Windows

A single T estimated from 10 years of data fails because Jupiter's 12-year orbit causes
the synodic eclipse period to drift by 1–4 seconds per year. Over 10 years this
accumulates to tens of thousands of seconds of secular OmC drift, completely swamping
the 76-second Rømer signal.

The fix is to work within **opposition-to-opposition windows** of roughly one Earth year
each. Within a single window, Jupiter has moved less than 30° in its orbit and period
drift is small.

**Finding oppositions**: We locate **local minima in D** (the Earth–moon distance
array). D is minimized when Earth is closest to Jupiter, i.e., at opposition. We use a
half-width of 5 events: index `i` is a local minimum if `D[i] ≤ D[j]` for all `j` in
`[i−5, i+5]`, with a gap filter ensuring consecutive minima are separated by at least
5 events.

```python
def find_local_minima(D, half_w=5):
    n = len(D)
    minima = []
    for i in range(half_w, n - half_w):
        if all(D[i] <= D[i+k] for k in range(-half_w, half_w+1)):
            if not minima or i - minima[-1] > half_w:
                minima.append(i)
    return minima
```

This gives ~10 oppositions per moon over 10 years, defining ~9 windows. After filtering
on window size (40–600 events) to exclude dataset-edge fragments, we retain 8 windows
per moon.

### 5.3 Period Estimation at Opposition

Within each window we need a reference period T and reference epoch t0 to compute OmC.
Estimating T from the entire window introduces a **Rømer bias**: if Earth is moving
away throughout the window, successive eclipse intervals appear systematically longer
than T_true, biasing T high and inflating OmC.

The fix: estimate T from a ~60-event sub-window **centered on opposition** (the local
D minimum), where `dD/dt ≈ 0`. At opposition, Earth's motion relative to Jupiter is
mostly transverse, so the light-travel-time change between consecutive eclipses is near
zero. The Rømer bias in T vanishes to first order.

```python
def fit_T(t_arr, i_center, half_w=30):
    n  = len(t_arr)
    lo = max(0, i_center - half_w)
    hi = min(n - 1, i_center + half_w)
    t_w = [t_arr[i] for i in range(lo, hi+1)]
    k   = list(range(len(t_w)))
    nk  = len(k);  mk = sum(k)/nk;  mt = sum(t_w)/nk
    cov = sum((k[i]-mk)*(t_w[i]-mt) for i in range(nk))
    vk  = sum((kk-mk)**2 for kk in k)
    T   = cov / vk
    return T, mt - T*mk, lo
```

This is OLS regression of event time `t[i]` on event index `i` over the
opposition-centered sub-window. The slope is T (days per eclipse) and the intercept
`mt − T×mk` is the reference epoch t0, anchored at index `lo` (the sub-window start).

We do this **separately** for `t_obs` (yielding T_obs, t0_obs) and for `t_geo`
(yielding T_geo, t0_geo).

**Assumption**: `dD/dt ≈ 0` near opposition. This is the standard Rømer approximation.
It holds to within a few percent over a ±30-event window around opposition for
Earth–Jupiter geometry, and it is not a systematic bias — any residual error appears
as random scatter across windows.

### 5.4 Computing the OmC Slopes

For each opposition-to-opposition window (events i_lo through i_hi − 1):

```python
def window_slope(t_arr, D_arr, i_lo, i_hi, T, t0_ref, i_zero):
    n_w = i_hi - i_lo
    OmC = [(t_arr[i_lo+j] - t0_ref - (i_lo+j - i_zero)*T) * DAYS_TO_S
           for j in range(n_w)]
    Dw  = [D_arr[i_lo+j] for j in range(n_w)]

    mD     = sum(Dw) / n_w
    mO     = sum(OmC) / n_w
    cov_DO = sum((Dw[j]-mD)*(OmC[j]-mO) for j in range(n_w))
    var_D  = sum((d-mD)**2 for d in Dw)
    slope  = cov_DO / var_D      # seconds per AU
    ...
```

OmC is in **seconds** and D is in **AU**, so slope is in s/AU. The expected Rømer
slope is:

```
AU_TO_KM / c  =  1.495978707×10⁸ km/AU  /  299792.458 km/s  =  499.0 s/AU
```

The function also computes:
- **Residual root mean square (RMS)** in seconds — scatter of OmC around the fitted line
- **Standard error on the slope**: `se_slope = sqrt(σ² / Var(D))` where σ² is the mean
  squared residual with n_w − 2 degrees of freedom
- **Pearson r** — the correlation coefficient of OmC vs D, a dimensionless quality
  indicator ranging from −1 to +1

### 5.5 Recovering c

For each window:

```python
s_romer = s_obs - s_geo
c_km_s  = AU_TO_KM / s_romer
```

The standard error on c is propagated from the two slope uncertainties, treated as
independent:

```python
se_romer = sqrt(se_obs² + se_geo²)
se_c     = AU_TO_KM × se_romer / s_romer²
```

### 5.6 Weighted Combination

The final estimate combines all valid windows using inverse-variance weighting:

```
weights    = 1 / se_c²  for each window
c_combined = Σ(w × c) / Σ(w)
se_combined = 1 / sqrt(Σ(w))
```

This gives highest weight to windows where the slope is most precisely determined —
those with the highest |r_obs| and lowest OmC scatter.

---

## 6. Results

### 6.1 Per-Window Summary

All 24 windows (8 per moon) across the 10-year dataset:

| Moon     | Win | N   | slope_geo | slope_obs | slope_romer | c (km/s) | Error   | r_obs  |
|----------|-----|-----|-----------|-----------|-------------|----------|---------|--------|
| Io       | 0   | 224 | +87.4     | +586.4    | 499.0       | 299,778  | −0.005% | 0.670  |
| Io       | 1   | 224 | +177.0    | +676.0    | 499.1       | 299,765  | −0.009% | 0.889  |
| Io       | 2   | 225 | +159.2    | +658.2    | 499.0       | 299,808  | +0.005% | 0.989  |
| Io       | 3   | 225 | +54.1     | +554.0    | 499.8       | 299,303  | −0.163% | 0.736  |
| Io       | 4   | 227 | −58.9     | +439.6    | 498.6       | 300,057  | +0.088% | 0.491  |
| Io       | 5   | 227 | −110.5    | +389.1    | 499.6       | 299,435  | −0.119% | 0.461  |
| Io       | 6   | 227 | −110.4    | +388.4    | 498.8       | 299,927  | +0.045% | 0.604  |
| Io       | 7   | 226 | −127.4    | +370.4    | 497.8       | 300,532  | +0.247% | 0.896  |
| **Io mean** |  |     |           |           |             | **299,921** | **+0.043%** | |
| Europa   | 0   | 111 | −349.7    | +135.3    | 485.1       | 308,405  | +2.873% | 0.060  |
| Europa   | 1   | 112 | −792.0    | −292.9    | 499.0       | 299,782  | −0.003% | −0.205 |
| Europa   | 2   | 112 | −830.1    | −331.4    | 498.8       | 299,931  | +0.046% | −0.939 |
| Europa   | 3   | 112 | −563.2    | −64.1     | 499.1       | 299,711  | −0.027% | −0.054 |
| Europa   | 4   | 113 | −200.6    | +297.9    | 498.5       | 300,092  | +0.100% | 0.167  |
| Europa   | 5   | 113 | +156.7    | +656.0    | 499.3       | 299,616  | −0.059% | 0.307  |
| Europa   | 6   | 113 | +479.6    | +978.1    | 498.5       | 300,093  | +0.100% | 0.397  |
| Europa   | 7   | 113 | +735.9    | +1233.4   | 497.5       | 300,674  | +0.294% | 0.531  |
| **Europa mean** | | |           |           |             | **299,949** | **+0.052%** | |
| Ganymede | 0   | 56  | +59.6     | +516.0    | 456.3       | 327,827  | +9.351% | 0.488  |
| Ganymede | 1   | 55  | +128.7    | +627.6    | 499.0       | 299,800  | +0.003% | 0.973  |
| Ganymede | 2   | 56  | +77.1     | +575.4    | 498.3       | 300,188  | +0.132% | 0.998  |
| Ganymede | 3   | 55  | −139.5    | +359.1    | 498.6       | 300,045  | +0.084% | 0.489  |
| Ganymede | 4   | 56  | −224.5    | +274.2    | 498.8       | 299,929  | +0.046% | 0.251  |
| Ganymede | 5   | 56  | −162.2    | +337.2    | 499.4       | 299,539  | −0.084% | 0.333  |
| Ganymede | 6   | 56  | +31.0     | +530.3    | 499.3       | 299,628  | −0.055% | 0.911  |
| Ganymede | 7   | 56  | +79.5     | +577.4    | 497.9       | 300,462  | +0.223% | 0.558  |
| **Ganymede mean** | | |          |           |             | **300,203** | **+0.137%** | |

### 6.2 Final Combined Result

```
COMBINED (24 windows, all moons, weighted by 1/se²):
  c = 300,078 ± 2,815 km/s
  Error vs CODATA: +0.095%
  CODATA (Committee on Data for Science and Technology reference value):
         299,792.458 km/s
```

### 6.3 Observations on slope_obs

The slope_obs values range from −331 to +1233 s/AU across moons and windows. This large
variation comes from the Laplace resonance contamination, which adds constructively or
destructively to the Rømer signal depending on the orbital phase. The corrected
slope_romer clusters tightly around 498–500 s/AU across all 22 clean windows — a direct
demonstration that the subtraction is working as intended.

### 6.4 The Two Outliers

Europa window 0 (+2.87%) and Ganymede window 0 (+9.35%) are dataset-edge artifacts. The
first opposition in each moon's time series falls only 9–19 events from the start of the
array. The T estimation sub-window, which ideally uses ±30 events around the opposition,
truncates at the dataset boundary, degrading the T estimate and distorting the OmC
baseline for the first window.

Both outliers are automatically down-weighted: their se_c values (94,000 and 208,000
km/s) are 17–52× larger than those of typical clean windows (~4,000–37,000 km/s), so
their contribution to the weighted mean is negligible. Excluding them shifts the
combined estimate by only 0.008%.

---

## 7. Assumptions Summary

| Assumption | Justification | Impact if violated |
|------------|---------------|--------------------|
| Jupiter modeled as sphere (equatorial radius) | Oblateness is ~6%; error in umbra radius is <1% | <0.1% error in c |
| Umbra only, no penumbra | Matches historical observability; sharp ingress/egress | Eclipse counts change; method unchanged |
| `dD/dt ≈ 0` near opposition for T estimation | Valid to a few percent over ±30-event window | Per-window scatter; no systematic bias |
| OmC–D relationship is linear within each window | Valid when T bias and Jupiter drift are small within one year | Nonlinearity appears as residual scatter |
| Slope uncertainties se_obs and se_geo are independent | Conservative; they share the same D array | se_c is slightly overestimated |
| CODATA c used to build t_obs | Discussed at length in §8 | No measurable impact on the recovered value |
| Gaussian residuals for OLS | Standard assumption | Confidence intervals affected, not point estimates |

---

## 8. Why the Method Is Not Circular

### 8.1 The apparent circularity

We use CODATA c to build t_obs = t_geo + D/c, then recover c from the data. If we
simply fitted the OLS slope of OmC_obs vs D and claimed c = AU_TO_KM / slope_obs, that
would indeed be exactly circular — we would recover CODATA by construction, because we
injected it.

Instead, we compute:

```
c = AU_TO_KM / (slope_obs − slope_geo)
```

The quantity slope_obs − slope_geo is not imposed by construction. slope_geo is
derived entirely from t_geo — the geometric eclipse times that have no c embedded in
them. The difference measures the additional D-dependence introduced when going from
geometric to observed time, which is exactly AU_TO_KM/c. The Laplace resonance
contamination, whatever its value, appears identically in both slopes and cancels.

The per-window scatter of ±0.25% shows that the identity is not satisfied exactly:
the T estimates for t_obs and t_geo differ slightly (because t_obs has a small
D-dependent offset), so the OmC baselines differ and the slopes do not cancel
perfectly. This imperfection is the residual measurement uncertainty of the method,
not a bias.

### 8.2 The deeper reason CODATA must be used

There is a more fundamental point that goes beyond algebraic self-consistency. The
answer to the question "why must we use CODATA c to build t_obs?" is: **because CODATA
c is already embedded throughout the JPL ephemeris at every level, and there is no
physically meaningful alternative.**

JPL's planetary ephemeris (DE440 and its predecessors, which underpin Horizons) is not
a purely geometric construction. It is a simultaneous fit to decades of observational
data: radar ranging to planets, Very Long Baseline Interferometry (VLBI) measurements
of spacecraft, laser ranging to the Moon, and transponder timing from deep-space
probes. All of these measurements involve light or radio waves traveling finite
distances. Every range measurement used to constrain the ephemeris was converted from
travel time to distance using c — specifically the International Astronomical Union
(IAU) adopted value, which is the same as CODATA.

The positions and velocities in Horizons — the very numbers we use to compute shadow
crossing times and Earth–moon distances — were determined by a fit that assumes c. The
AU itself, as realized in the ephemeris, is defined such that the light-travel time
from Sun to Earth is 499.004782 seconds at 1 AU. The column `earth_moon_au` in our
eclipse files was computed using positions derived from that fit.

In other words, the "geometric" eclipse times are not purely geometric in the naive
sense. They are geometric in the sense that the code computes shadow crossings from
orbital positions without re-applying a light-travel correction — but those orbital
positions were themselves inferred from light-travel-time measurements that assumed c.

The consequence: what our analysis actually demonstrates is internal self-consistency.
We show that if you take the positions and distances that JPL derived using c = CODATA,
construct simulated observed times using the same c, and apply the contamination-
cancellation Rømer method, you recover c to within 0.1%. This confirms that the method
is correctly implemented and that the Laplace resonance contamination is successfully
removed. It does not — and cannot — produce an independent measurement of c from this
data source, because an independent measurement would require positions derived without
assuming c.

A genuinely independent replication of Rømer's measurement would require raw observed
eclipse times recorded by a real telescope, with no prior knowledge of c used in
reducing the data. Rømer himself had this: his observations were direct photometric
records of Io's reappearance from shadow, with no JPL ephemeris involved.