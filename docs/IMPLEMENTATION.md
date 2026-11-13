# Implementation Details

This document explains how the simulation code implements the system model and the algorithms of the paper. Equation numbers (1)–(28) refer to the manuscript.

---

## 1. Overview of the simulation loop

Both simulators (`src/uet_ris_sim.py` for one vehicle, `src/mu_sim.py` for several) run the same frame-by-frame loop. Frames are 10 ms long. For every frame *k* and every scheme:

```
1. EKF prediction            s_pred, P_pred                         (10)
2. Decision                  keep current RIS phases, or reconfigure  (scheme-specific)
3. Communication             SNR at the TRUE vehicle position         (4) / (23)
                             effective rate with overhead charged     (6) / (25)
4. Sensing                   echo measurement whose noise depends on
                             the illumination of the applied beam     (8), (9)
5. EKF update                s_est, P_est                             (11)
```

All schemes in one trial share the same vehicle trajectory, initial EKF state and RIS geometry, so the comparison is paired. A trial returns the time-averaged rate, outage, reconfiguration frequency, tracking RMSE and decision time of every scheme. `aggregate` / `agg` then average over trials and compute 95% confidence intervals.

---

## 2. Geometry and channel (Sec. 3.1–3.2)

| Item | Code | Notes |
|---|---|---|
| RSU position, RIS centre | `P_B = (-50, 40, 12)`, `P_R = (0, 0, 6)` | metres; the RIS lies in the wall plane y = 0 and faces +y |
| Element positions, eq. (1) | `ris_elements(nx, nz)` | half-wavelength UPA in the x–z plane; also returns each element's column index |
| RSU–RIS coefficients *g_n*, eq. (2) | `RIS.__init__` → `amp1`, `ph1` | includes √(P_t N_t) (MRT array gain) and the element gain *G_e* |
| RIS–vehicle coefficients *h_n(p)* and cascade *c_n(p)*, eq. (2), (4) | `RIS.coeffs(P)` | **exact spherical wavefront**: distance from every element to the 3-D vehicle antenna (height `Z_V` = 1.5 m), vectorised over many positions |
| SNR γ(θ, p), eq. (4) | `RIS.snr(phi, P)` | `|c @ exp(j·phi)|² / NOISE_W` |
| Continuous-phase maximum γ^max | `RIS.snr_max(P)` | `(Σ|c_n|)² / NOISE_W`, used for the illumination ratio η |
| Reconfiguration overhead τ_c, eq. (5) | `RIS.__init__` → `tau_c = N*BITS/R_CTRL + T_SW` | ρ = τ_c / T is `ovh` in `run_trial` |

The direct RSU–vehicle path is assumed blocked, so only the RIS path is simulated.

---

## 3. Phase design (Sec. 5.2)

`RIS.focus(Q)` implements eq. (16):

* `Q` holds one focal point (ordinary focusing) or two (the **segmented** candidate).
* With two points, the RIS columns are split into contiguous halves. Each half compensates the full RSU→element→focal-point phase for its own point.
* Phases are wrapped to [0, 2π) and quantised to the *b*-bit alphabet *F_b* (`BITS = 2` → 4 levels).

The multi-vehicle version `MURIS.focus_sub(phi, k, Q)` does the same, but only for the elements of sub-array *k*, and leaves every other element unchanged.

---

## 4. Tracking (Sec. 3.4)

| Item | Code |
|---|---|
| State **s** = [x, y, vx, vy], constant-velocity model, eq. (7) | `F_mat(T)`, `Q_mat(T, SIG_A)` |
| Measurement **h(s)** = [range from RIS, azimuth from RIS normal], eq. (8) | `h_meas(s)` |
| Jacobian **H** | `H_jac(s)` |
| Sensing–illumination coupling, eq. (9) | `sens_std(d, g)`: σ scales as (d/D0)² / √η, floored at η = 10⁻³ and capped at 5 m / 10° |
| EKF predict/update, eqs. (10)–(11) | inline in `run_trial` / `run_mu`; the angle innovation is wrapped to (−π, π] |

The measurement noise actually drawn uses η computed from the **applied** beam at the **true** position. A misaligned beam therefore produces noisier echoes, which is the coupling discussed in the paper.

Trajectories (`trajectory`, `mu_traj`) start near x = −45 m (or at a random position and direction in the multi-vehicle case). Speeds are drawn within ±10% of nominal with random longitudinal acceleration (std 1 m/s²), and half of the single-vehicle trajectories include a smooth 3.5 m lane change lasting 2 s.

---

## 5. Single-vehicle schemes (method names in `run_trial`)

| Method | What it does |
|---|---|
| `genie` | Focuses on the true position in every frame with no overhead (upper reference) |
| `periodicM` (M = 1, 5, 10, 20) | Every M frames, focuses on the predicted position at the centre of the hold window, (M−1)/2 frames ahead; overhead charged in those frames |
| `et_mean` | Proposed trigger evaluated only at the predicted mean (ablation) |
| `uet` | **Proposed UET-RIS** (below) |
| `random` | Fixed random phases (lower reference) |

### UET-RIS step by step (Algorithm 1 / Table 3)

1. **Sigma points**, eq. (14): `sigma_points(m, P)` returns 5 points with weights κ/(n+κ) = 1/3 and 1/6 (n = 2, κ = 1).
2. **Expected rate**, eq. (15): `exp_rate(ph)` = Σ wᵢ log₂(1 + γ(ph, χᵢ)).
3. **Fresh candidates**: focusing at the predicted mean, plus the segmented candidate at mean ± 0.8 σ∥ along the velocity direction. `Rf` is the best expected rate and `Ro` that of the current phases.
4. **Regret accumulation**, eq. (17): `acc += max(Rf - Ro, 0)`.
5. **Trigger**, eq. (18): reconfigure if `acc >= ovh * Rf` (ski-rental break-even), **or** if the worst sigma point of the current phases is below `GAMMA_MIN_DB` = 3 dB and a fresh profile is better (the reliability guard).
6. **Lead-compensated design**, eq. (19): lead = round(H̄/2), clipped to `LEAD_MAX` = 5. The state and covariance are propagated by F^lead, both candidates are rebuilt at the lead point, and the better one is applied. Then `acc = 0` and H̄ ← ½H̄ + ½·(frames since the last update).

The decision time of every scheme is measured with `time.perf_counter()` around step 2 of the loop and reported as `tdec_ms`.

---

## 6. Multi-vehicle extension (Sec. 6, `src/mu_sim.py`)

| Item | Code |
|---|---|
| *K* column sub-arrays | `MURIS` → `grp`, `masks` (contiguous groups of `nx/K` columns) |
| Per-vehicle SNR over the whole surface, eq. (23) | `ris.coeffs(p_k) @ exp(j·phi)` with the full phase vector |
| Event cost, eq. (24) | `ovh = (len(upd) * nbits_sub / R_CTRL + tsw) / T_FRAME` if any sub-array is updated |
| Rate with shared overhead, eq. (25) | `rates[t, k] = (1 - ovh) * log2(1 + g)` |
| Set-up and item costs, eq. (26) | `c_setup = rho_0 * Rsum`, `c_item = rho_m * Rsum` |
| Event trigger, eq. (27) | `(A >= c_setup + c_item).any() or guard.any()` |
| Batch, eq. (28) | `{argmax A} ∪ {k : A_k >= XI * c_item} ∪ {guarded}`, with `XI = 0.5` |
| Independent baseline `indep` | each vehicle updates when its own `A_k >= c_setup + c_item` |

Sensing for vehicle *k* uses the illumination ratio relative to its **own sub-array's** maximum (`snrmax_sub`). Each trial lasts `T_SIM` = 200 frames.

---

## 7. Experiments and seeds

| Experiment | Script call | Sweep | Trials/point | Seeds |
|---|---|---|---|---|
| Rate vs. speed (Figs. 3, 5) | `uet_ris_sim.py` | v = 10…40 m/s, N = 1024 | 40 | `1000·v + i` |
| Rate vs. RIS size (Fig. 4) | `uet_ris_sim.py` | N = 64…1024, v = 30 m/s | 40 | `7·nx + 100000 + i` |
| Sensing accuracy (Fig. 6) | `uet_ris_sim.py` | σθ,0 = 0.5…4°, v = 25 m/s | 40 | `500000 + 1000·σ + i` |
| Trace (Fig. 7) | `uet_ris_sim.py` | v = 30 m/s | 1 | `0` |
| Hold-time law (Fig. 8) | `scaling.py` | post-processes `results.json` | – | – |
| Multi-vehicle vs. K (Fig. 9) | `mu_sim.py` | K = 2, 4, 8 | 30 | `7000 + 100·K + i` |
| Multi-vehicle vs. speed (Table 6) | `mu_sim.py` | v = 10…40 m/s, K = 4 | 30 | `9000 + 100·v + i` |
| Multi-vehicle vs. T_sw (Fig. 10) | `mu_sim.py` | T_sw = 0.1…1 ms, K = 4 | 30 | `20000 + i` (common across T_sw) |

Trials run in parallel with `multiprocessing.Pool` on all available CPU cores.

---

## 8. Output files

`results/results.json`
```
{ "speed": {v:  {method: {rate, outage, recfg, rmse, tdec_ms, snr_db, rate_ci95}}},
  "size":  {nx: {...}},  "sens": {sigma_deg: {...}},
  "trace": {method: {trace_rate, trace_err, trace_snr}, "traj": [[x, y], ...]},
  "MC": 40, "params": {...} }
```

`results/mu_results.json`
```
{ "K": {K: {method: {rate, jain, events, bits, minrate, outage, tdec_ms, rate_ci95}}},
  "v": {...}, "tsw": {T_sw_ms: {...}}, "MC": 30, "XI": 0.5, "T_SIM": 200 }
```

Field meanings: `rate` = average effective spectral efficiency (bit/s/Hz); `recfg` / `events` = reconfiguration (switching) events per frame; `bits` = control bits per frame; `outage` = fraction of frames with SNR < 3 dB; `jain` = Jain's fairness index over vehicles; `rmse` = EKF position RMSE (m); `tdec_ms` = mean decision time per frame (ms).

---

## 9. Changing parameters and extending the code

* Global physical parameters are module-level constants at the top of `src/uet_ris_sim.py` (see README Sec. 5). `mu_sim.py` imports them, so a single change affects both simulators.
* To add a new single-vehicle scheme, add its name to `METHODS` and a branch to the decision block of `run_trial` that sets `phi` and `recfg`.
* To use a different RIS size in the multi-vehicle study, change `MURIS(32, 32, K)` in `run_mu`. `nx` must be divisible by `K`.
* The quick reproduction (`python reproduce.py --quick`) is a good regression test after any change.

---

## 10. Known limitations

* The direct RSU–vehicle path is ignored (assumed blocked), and the RSU–RIS link is treated as a static LoS channel.
* Sensing accuracy follows the inverse-SNR scaling of eq. (9) rather than a full waveform-level radar simulation.
* The decision times are from unoptimised Python. A compiled implementation would be needed for K = 8 on embedded hardware.
* In the multi-vehicle case the "genie" scheme is a reference, not an upper bound (see Sec. 8.7 of the paper).
