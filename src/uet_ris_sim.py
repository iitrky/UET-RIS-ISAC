"""
Simulation for: Uncertainty-Aware Event-Triggered Reconfiguration (UET-RIS)
for semi-passive RIS-aided ISAC in real-time vehicular (V2I) links.

Everything reported in the manuscript is produced by this script.
Run:  python3 uet_ris_sim.py   (writes results/*.json and figures/*.png)
"""
import numpy as np, json, time, os, sys
from multiprocessing import Pool

# ----------------------------------------------------------------- parameters
C = 3e8; FC = 28e9; LAM = C / FC; K0 = 2 * np.pi / LAM
BW = 50e6; NF_DB = 7.0
NOISE_W = 10 ** ((-174 + 10 * np.log10(BW) + NF_DB - 30) / 10)
PT_W = 10 ** ((33 - 30) / 10)          # 33 dBm
NT = 32                                # RSU antennas (MRT towards RIS)
GE = 4.0                               # RIS element gain (6 dB, product of both sides)
T_FRAME = 10e-3                        # frame length [s]
R_CTRL = 1e6                           # RIS control-link rate [bit/s]
T_SW = 0.2e-3                          # RIS switching/settling time [s]
BITS = 2                               # phase-shifter resolution
GAMMA_MIN_DB = 3.0                     # reliability guard (SNR threshold)
LEAD_MAX = 5                           # max design lead (frames)
SIG_A = 1.0                            # EKF process noise (accel std) [m/s^2]
SIG_R0 = 0.10; SIG_TH0 = np.deg2rad(1.0); D0 = 20.0   # sensing accuracy at D0, full gain
P_B = np.array([-50.0, 40.0, 12.0])    # RSU position
P_R = np.array([0.0, 0.0, 6.0])        # RIS centre (wall y=0, normal +y)
Z_V = 1.5                              # vehicle antenna height
X_START, X_END = -45.0, 45.0


def ris_elements(nx, nz):
    d = LAM / 2
    xs = (np.arange(nx) - (nx - 1) / 2) * d
    zs = (np.arange(nz) - (nz - 1) / 2) * d
    X, Z = np.meshgrid(xs, zs, indexing="ij")          # column index = x
    pos = np.stack([X.ravel(), np.zeros(X.size), Z.ravel()], 1) + P_R
    col = np.repeat(np.arange(nx), nz)
    return pos, col


class RIS:
    def __init__(self, nx, nz):
        self.pos, self.col = ris_elements(nx, nz)
        self.nx, self.N = nx, nx * nz
        self.d1 = np.linalg.norm(self.pos - P_B, axis=1)
        self.amp1 = np.sqrt(PT_W * NT) * GE * (LAM / (4 * np.pi)) ** 2 / self.d1
        self.ph1 = K0 * self.d1
        self.tau_c = self.N * BITS / R_CTRL + T_SW       # reconfiguration overhead [s]

    def coeffs(self, P):
        """cascaded per-element coefficients for vehicle positions P (M,2)->(M,N)"""
        P3 = np.column_stack([P, np.full(len(P), Z_V)])
        d2 = np.linalg.norm(P3[:, None, :] - self.pos[None], axis=2)
        return (self.amp1 / d2) * np.exp(-1j * (self.ph1 + K0 * d2))

    def snr(self, phi, P):
        c = self.coeffs(np.atleast_2d(P))
        return np.abs(c @ np.exp(1j * phi)) ** 2 / NOISE_W

    def snr_max(self, P):
        c = self.coeffs(np.atleast_2d(P))
        return np.sum(np.abs(c), 1) ** 2 / NOISE_W

    def focus(self, Q):
        """segmented focusing: Q (K,2) focal points, contiguous column segments"""
        K = len(Q)
        seg = np.minimum((self.col * K) // self.nx, K - 1)
        Q3 = np.column_stack([Q, np.full(K, Z_V)])
        d2 = np.linalg.norm(self.pos - Q3[seg], axis=1)
        phi = np.mod(self.ph1 + K0 * d2, 2 * np.pi)
        step = 2 * np.pi / 2 ** BITS
        return np.mod(np.round(phi / step) * step, 2 * np.pi)


# ----------------------------------------------------------------- tracking
def F_mat(T):
    return np.array([[1, 0, T, 0], [0, 1, 0, T], [0, 0, 1, 0], [0, 0, 0, 1.0]])


def Q_mat(T, s):
    q = np.array([[T**4 / 4, T**3 / 2], [T**3 / 2, T**2]]) * s**2
    Q = np.zeros((4, 4)); Q[np.ix_([0, 2], [0, 2])] = q; Q[np.ix_([1, 3], [1, 3])] = q
    return Q


def h_meas(s):
    dx, dy, dz = s[0] - P_R[0], s[1] - P_R[1], Z_V - P_R[2]
    r = np.sqrt(dx**2 + dy**2 + dz**2)
    return np.array([r, np.arctan2(dx, dy)])


def H_jac(s):
    dx, dy, dz = s[0] - P_R[0], s[1] - P_R[1], Z_V - P_R[2]
    r = np.sqrt(dx**2 + dy**2 + dz**2); rho2 = dx**2 + dy**2
    return np.array([[dx / r, dy / r, 0, 0], [dy / rho2, -dx / rho2, 0, 0]])


def sens_std(d, g):
    g = max(g, 1e-3); sc = (d / D0) ** 2 / np.sqrt(g)
    return min(SIG_R0 * sc, 5.0), min(SIG_TH0 * sc, np.deg2rad(10))


def sigma_points(m, P):
    n, kap = 2, 1.0
    S = np.linalg.cholesky((n + kap) * (P + 1e-9 * np.eye(2)))
    pts = np.vstack([m, m + S.T, m - S.T])
    w = np.array([kap / (n + kap)] + [1 / (2 * (n + kap))] * 4)
    return pts, w


# ----------------------------------------------------------------- trajectory
def trajectory(v_nom, rng):
    v = v_nom * rng.uniform(0.9, 1.1)
    y0 = rng.choice([12.0, 15.5]); lc = rng.random() < 0.5
    x, vx, xs, ys = X_START + rng.uniform(-3, 3), v, [], []
    t, t_lc = 0.0, rng.uniform(0.5, 1.5) * (X_END - X_START) / v / 2
    while x < X_END:
        vx = max(3.0, vx + rng.normal(0, 1.0) * T_FRAME)
        x += vx * T_FRAME; t += T_FRAME
        y = y0
        if lc:   # smooth 3.5 m lane change over 2 s
            u = np.clip((t - t_lc) / 2.0, 0, 1)
            y = y0 + (3.5 if y0 == 12.0 else -3.5) * (u - np.sin(2 * np.pi * u) / (2 * np.pi))
        xs.append(x); ys.append(y)
    return np.column_stack([xs, ys])


# ----------------------------------------------------------------- one run
METHODS = ["genie", "periodic1", "periodic5", "periodic10", "periodic20", "et_mean", "uet", "random"]


def run_trial(args):
    global SIG_TH0, SIG_R0, R_CTRL
    v_nom, nx, seed = args[:3]
    th = args[3] if len(args) > 3 else 1.0
    R_CTRL = args[4] if len(args) > 4 else 1e6
    SIG_TH0, SIG_R0 = np.deg2rad(th), 0.1 * th
    rng = np.random.default_rng(seed)
    ris = RIS(nx, nx)
    traj = trajectory(v_nom, rng)
    K_ = len(traj)
    F, Q = F_mat(T_FRAME), Q_mat(T_FRAME, SIG_A)
    vel0 = (traj[1] - traj[0]) / T_FRAME
    s0 = np.r_[traj[0] + rng.normal(0, 1, 2), vel0 + rng.normal(0, 1, 2)]
    P0 = np.diag([1, 1, 1, 1.0])
    gmin = 10 ** (GAMMA_MIN_DB / 10)
    ovh = ris.tau_c / T_FRAME
    snrmax = ris.snr_max(traj)
    phi_rand = rng.uniform(0, 2 * np.pi, ris.N)
    out = {}
    for m in METHODS:
        s, P = s0.copy(), P0.copy()
        phi = ris.focus(s[None, :2]); nre = 0; acc = 0.0; hold = 4.0; since = 0
        rate, snr_l, err, tdec = np.zeros(K_), np.zeros(K_), np.zeros(K_), []
        for k in range(K_):
            sp = F @ s; Pp = F @ P @ F.T + Q
            t0 = time.perf_counter(); recfg = False
            if m == "genie":
                phi = ris.focus(traj[k][None]); recfg = False   # zero overhead bound
            elif m.startswith("periodic"):
                M = int(m[8:])
                if k % M == 0:      # design for the centre of the hold window
                    sc = np.linalg.matrix_power(F, (M - 1) // 2) @ sp
                    phi = ris.focus(sc[None, :2]); recfg = True
            elif m in ("et_mean", "uet"):
                # --- uncertainty-aware evaluation set (sigma points or mean only)
                if m == "uet":
                    pts, W = sigma_points(sp[:2], Pp[:2, :2])
                else:
                    pts, W = sp[None, :2], np.ones(1)
                c = ris.coeffs(pts)
                def exp_rate(ph):
                    g = np.abs(c @ np.exp(1j * ph)) ** 2 / NOISE_W
                    return float(np.log2(1 + g) @ W), g
                cands = [sp[None, :2]]
                if m == "uet":   # segmented (beam-broadened) candidate along heading
                    vd = sp[2:] / max(np.linalg.norm(sp[2:]), 1e-6)
                    sa = np.sqrt(max(vd @ Pp[:2, :2] @ vd, 1e-8))
                    cands.append(np.array([sp[:2] - 0.8 * sa * vd, sp[:2] + 0.8 * sa * vd]))
                Rf = max(exp_rate(ris.focus(q))[0] for q in cands)
                Ro, go = exp_rate(phi)
                acc += max(Rf - Ro, 0.0)                  # accumulated regret (rent)
                guard = (go.min() < gmin) and (Rf > Ro)
                if acc >= ovh * Rf or guard:               # ski-rental break-even
                    lead = int(np.clip(round(hold / 2), 0, LEAD_MAX))
                    Fl = np.linalg.matrix_power(F, lead)
                    sl = Fl @ sp; Pl = Fl @ Pp @ Fl.T
                    if m == "uet":
                        pts_l, W_l = sigma_points(sl[:2], Pl[:2, :2])
                        vd = sl[2:] / max(np.linalg.norm(sl[2:]), 1e-6)
                        sa = np.sqrt(max(vd @ Pl[:2, :2] @ vd, 1e-8))
                        cl = [sl[None, :2], np.array([sl[:2] - 0.8 * sa * vd, sl[:2] + 0.8 * sa * vd])]
                    else:
                        pts_l, W_l = sl[None, :2], np.ones(1); cl = [sl[None, :2]]
                    c_l = ris.coeffs(pts_l); bestv = -1
                    for q in cl:
                        ph = ris.focus(q)
                        v_ = float(np.log2(1 + np.abs(c_l @ np.exp(1j * ph)) ** 2 / NOISE_W) @ W_l)
                        if v_ > bestv: bestv, phi_new = v_, ph
                    phi = phi_new; recfg = True; acc = 0.0
                    hold = 0.5 * hold + 0.5 * since; since = 0
                since += 1
            elif m == "random":
                phi = phi_rand
            tdec.append(time.perf_counter() - t0)
            nre += recfg
            g = ris.snr(phi, traj[k])[0]
            snr_l[k] = g
            rate[k] = (1 - ovh * recfg) * np.log2(1 + g)
            # ISAC echo measurement (accuracy depends on illumination gain)
            z_true = h_meas(np.r_[traj[k], 0, 0])
            sr, sth = sens_std(z_true[0], g / snrmax[k])
            z = z_true + rng.normal(0, 1, 2) * [sr, sth]
            Hk = H_jac(sp); Rk = np.diag([sr**2, sth**2])
            y = z - h_meas(sp); y[1] = (y[1] + np.pi) % (2 * np.pi) - np.pi
            S = Hk @ Pp @ Hk.T + Rk; Kg = Pp @ Hk.T @ np.linalg.inv(S)
            s = sp + Kg @ y; P = (np.eye(4) - Kg @ Hk) @ Pp
            err[k] = np.linalg.norm(s[:2] - traj[k])
        out[m] = dict(rate=rate.mean(), outage=float(np.mean(10*np.log10(snr_l) < GAMMA_MIN_DB)),
                      recfg=nre / K_, rmse=float(np.sqrt(np.mean(err**2))),
                      tdec_ms=1e3 * float(np.mean(tdec)), snr_db=float(np.mean(10*np.log10(snr_l))))
        if seed == 0:
            out[m]["trace_rate"] = rate.tolist(); out[m]["trace_err"] = err.tolist()
            out[m]["trace_snr"] = (10*np.log10(snr_l)).tolist()
    if seed == 0:
        out["traj"] = traj.tolist()
    return out


def aggregate(res):
    agg = {}
    for m in METHODS:
        agg[m] = {key: float(np.mean([r[m][key] for r in res]))
                  for key in ["rate", "outage", "recfg", "rmse", "tdec_ms", "snr_db"]}
        agg[m]["rate_ci95"] = float(1.96 * np.std([r[m]["rate"] for r in res]) / np.sqrt(len(res)))
    return agg


if __name__ == "__main__":
    MC = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    os.makedirs("results", exist_ok=True)
    speeds = [10, 15, 20, 25, 30, 35, 40]
    sizes = [8, 12, 16, 24, 32]
    jobs = {}
    for v in speeds:
        jobs[("speed", v)] = [(v, 32, 1000 * v + i) for i in range(MC)]
    for nx in sizes:
        jobs[("size", nx)] = [(30, nx, 7 * nx + 100000 + i) for i in range(MC)]
    sens = [0.5, 1.0, 2.0, 3.0, 4.0]
    for th in sens:
        jobs[("sens", th)] = [(25, 32, int(500000 + 1000 * th) + i, th) for i in range(MC)]
    jobs[("trace", 0)] = [(30, 32, 0)]
    flat = [(key, a) for key, lst in jobs.items() for a in lst]
    t0 = time.time()
    with Pool(os.cpu_count()) as pool:
        outs = pool.map(run_trial, [a for _, a in flat])
    print(f"simulated {len(flat)} trials in {time.time()-t0:.0f}s")
    grouped = {}
    for (key, _), o in zip(flat, outs):
        grouped.setdefault(key, []).append(o)
    final = {"speed": {v: aggregate(grouped[("speed", v)]) for v in speeds},
             "size": {nx: aggregate(grouped[("size", nx)]) for nx in sizes},
             "sens": {th: aggregate(grouped[("sens", th)]) for th in sens},
             "trace": grouped[("trace", 0)][0], "MC": MC,
             "params": dict(fc_GHz=FC/1e9, BW_MHz=BW/1e6, Pt_dBm=33, NT=NT, bits=BITS,
                            T_ms=T_FRAME*1e3, Rctrl_Mbps=R_CTRL/1e6, Tsw_ms=T_SW*1e3,
                            noise_dBm=10*np.log10(NOISE_W)+30)}
    with open("results/results.json", "w") as f:
        json.dump(final, f)
    for v in speeds:
        a = final["speed"][v]
        print(v, {m: (round(a[m]["rate"], 2), round(a[m]["recfg"], 3)) for m in METHODS})
    for nx in sizes:
        a = final["size"][nx]
        print(nx*nx, {m: (round(a[m]["rate"], 2), round(a[m]["recfg"], 3)) for m in METHODS})
