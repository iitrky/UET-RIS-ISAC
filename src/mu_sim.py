"""
Multi-vehicle extension: K vehicles share one semi-passive RIS (column sub-arrays, one per vehicle,
orthogonal sub-bands).  A reconfiguration *event* updates a set S of sub-arrays and costs
    tau(S) = |S| (N/K) b / R_c + T_sw        (shared switching/settling time T_sw)
of airtime for ALL vehicles in that frame.  This is an online joint-replenishment problem (JRP):
T_sw is the shared set-up cost, the per-sub-array bits are item costs.

Schemes: genie, periodic1, periodic5, indep (per-vehicle ski-rental, events merge only by chance),
         jrp (proposed: batched regret trigger with opportunistic inclusion).
Run: python3 mu_sim.py [MC]   ->  results/mu_results.json
"""
import numpy as np, json, os, sys, time
from multiprocessing import Pool
import uet_ris_sim as U

METHODS = ["genie", "periodic1", "periodic5", "indep", "jrp"]
XI = 0.5          # opportunistic inclusion threshold (fraction of a sub-array's item cost)
T_SIM = 200       # frames per trial (2 s)


class MURIS(U.RIS):
    def __init__(self, nx, nz, K):
        super().__init__(nx, nz)
        self.K = K
        self.grp = np.minimum((self.col * K) // nx, K - 1)          # sub-array of each element
        self.masks = [self.grp == k for k in range(K)]
        self.cols_k = nx // K

    def focus_sub(self, phi, k, Q):
        """return copy of phi with sub-array k focused on points Q (1 or 2 segments)"""
        m = self.masks[k]; out = phi.copy()
        pos = self.pos[m]; colk = self.col[m] - self.col[m].min()
        seg = np.minimum((colk * len(Q)) // self.cols_k, len(Q) - 1)
        Q3 = np.column_stack([Q, np.full(len(Q), U.Z_V)])
        d2 = np.linalg.norm(pos - Q3[seg], axis=1)
        ph = np.mod(self.ph1[m] + U.K0 * d2, 2 * np.pi)
        step = 2 * np.pi / 2 ** U.BITS
        out[m] = np.mod(np.round(ph / step) * step, 2 * np.pi)
        return out


def mu_traj(K, v_nom, rng, n):
    trajs = []
    for k in range(K):
        v = v_nom * rng.uniform(0.9, 1.1); d = rng.choice([-1, 1])
        y = rng.choice([12.0, 15.5])
        x = rng.uniform(-45, 0) if d > 0 else rng.uniform(0, 45)
        vx = d * v; pts = []
        for _ in range(n):
            vx = d * max(3.0, abs(vx) + rng.normal(0, 1.0) * U.T_FRAME)
            x += vx * U.T_FRAME; pts.append((x, y))
        trajs.append(np.array(pts))
    return trajs


def run_mu(args):
    K, v_nom, seed, tsw = args
    U.T_SW = tsw
    rng = np.random.default_rng(seed)
    ris = MURIS(32, 32, K)
    trajs = mu_traj(K, v_nom, rng, T_SIM)
    F, Q = U.F_mat(U.T_FRAME), U.Q_mat(U.T_FRAME, U.SIG_A)
    nbits_sub = ris.N // K * U.BITS
    rho_m = nbits_sub / U.R_CTRL / U.T_FRAME             # marginal (item) overhead fraction
    rho_0 = tsw / U.T_FRAME                              # shared set-up overhead fraction
    gmin = 10 ** (U.GAMMA_MIN_DB / 10)
    # sub-array-only maximum SNR (for sensing illumination ratio)
    def snrmax_sub(k, p):
        c = ris.coeffs(p[None])[0][ris.masks[k]]
        return np.sum(np.abs(c)) ** 2 / U.NOISE_W
    s0 = []
    for k in range(K):
        v0 = (trajs[k][1] - trajs[k][0]) / U.T_FRAME
        s0.append(np.r_[trajs[k][0] + rng.normal(0, 1, 2), v0 + rng.normal(0, 1, 2)])
    out = {}
    for m in METHODS:
        S = [x.copy() for x in s0]; P = [np.eye(4) for _ in range(K)]
        phi = np.zeros(ris.N)
        for k in range(K): phi = ris.focus_sub(phi, k, S[k][None, :2])
        A = np.zeros(K); hold = np.full(K, 4.0); since = np.zeros(K)
        rates = np.zeros((T_SIM, K)); snr_db = np.zeros((T_SIM, K))
        events = 0; bits = 0; tdec = []
        for t in range(T_SIM):
            sp = [F @ S[k] for k in range(K)]; Pp = [F @ P[k] @ F.T + Q for k in range(K)]
            t0 = time.perf_counter(); upd = []
            if m == "genie":
                for k in range(K): phi = ris.focus_sub(phi, k, trajs[k][t][None])
            elif m == "periodic1" or (m == "periodic5" and t % 5 == 0):
                M = 1 if m == "periodic1" else 5
                for k in range(K):
                    sc = np.linalg.matrix_power(F, (M - 1) // 2) @ sp[k]
                    phi = ris.focus_sub(phi, k, sc[None, :2])
                upd = list(range(K))
            elif m in ("indep", "jrp"):
                Rf, Ro, guard, cand = np.zeros(K), np.zeros(K), np.zeros(K, bool), [None] * K
                for k in range(K):
                    pts, W = U.sigma_points(sp[k][:2], Pp[k][:2, :2])
                    c = ris.coeffs(pts)
                    er = lambda ph: (np.log2(1 + np.abs(c @ np.exp(1j * ph)) ** 2 / U.NOISE_W) @ W,
                                     np.abs(c @ np.exp(1j * ph)) ** 2 / U.NOISE_W)
                    Ro[k], go = er(phi)
                    vd = sp[k][2:] / max(np.linalg.norm(sp[k][2:]), 1e-6)
                    sa = np.sqrt(max(vd @ Pp[k][:2, :2] @ vd, 1e-8))
                    Rf[k] = max(er(ris.focus_sub(phi, k, q))[0] for q in
                                (sp[k][None, :2], np.array([sp[k][:2] - .8 * sa * vd, sp[k][:2] + .8 * sa * vd])))
                    A[k] += max(Rf[k] - Ro[k], 0.0)
                    guard[k] = (go.min() < gmin) and (Rf[k] > Ro[k])
                Rsum = Rf.sum()
                c_item, c_setup = rho_m * Rsum, rho_0 * Rsum
                if m == "indep":
                    upd = [k for k in range(K) if A[k] >= c_setup + c_item or guard[k]]
                else:   # JRP-style batched trigger
                    if (A >= c_setup + c_item).any() or guard.any():
                        upd = sorted(set([k for k in range(K) if A[k] >= XI * c_item or guard[k]]
                                         + [int(np.argmax(A))]))
                for k in upd:   # lead-compensated design for each updated sub-array
                    lead = int(np.clip(round(hold[k] / 2), 0, U.LEAD_MAX))
                    Fl = np.linalg.matrix_power(F, lead); sl = Fl @ sp[k]; Pl = Fl @ Pp[k] @ Fl.T
                    pts, W = U.sigma_points(sl[:2], Pl[:2, :2]); c = ris.coeffs(pts)
                    vd = sl[2:] / max(np.linalg.norm(sl[2:]), 1e-6)
                    sa = np.sqrt(max(vd @ Pl[:2, :2] @ vd, 1e-8)); best, bv = None, -1
                    for q in (sl[None, :2], np.array([sl[:2] - .8 * sa * vd, sl[:2] + .8 * sa * vd])):
                        ph = ris.focus_sub(phi, k, q)
                        val = np.log2(1 + np.abs(c @ np.exp(1j * ph)) ** 2 / U.NOISE_W) @ W
                        if val > bv: best, bv = ph, val
                    phi = best; A[k] = 0.0; hold[k] = 0.5 * hold[k] + 0.5 * since[k]; since[k] = 0
                since += 1
            tdec.append(time.perf_counter() - t0)
            ovh = (len(upd) * nbits_sub / U.R_CTRL + tsw) / U.T_FRAME if upd else 0.0
            events += bool(upd); bits += len(upd) * nbits_sub
            for k in range(K):
                g = np.abs(ris.coeffs(trajs[k][t][None])[0] @ np.exp(1j * phi)) ** 2 / U.NOISE_W
                snr_db[t, k] = 10 * np.log10(g); rates[t, k] = (1 - ovh) * np.log2(1 + g)
                z_true = U.h_meas(np.r_[trajs[k][t], 0, 0])
                sr, sth = U.sens_std(z_true[0], g / snrmax_sub(k, trajs[k][t]))
                z = z_true + rng.normal(0, 1, 2) * [sr, sth]
                Hk = U.H_jac(sp[k]); Rk = np.diag([sr**2, sth**2])
                y = z - U.h_meas(sp[k]); y[1] = (y[1] + np.pi) % (2 * np.pi) - np.pi
                Sk = Hk @ Pp[k] @ Hk.T + Rk; Kg = Pp[k] @ Hk.T @ np.linalg.inv(Sk)
                S[k] = sp[k] + Kg @ y; P[k] = (np.eye(4) - Kg @ Hk) @ Pp[k]
        ur = rates.mean(0)
        out[m] = dict(rate=float(ur.mean()), jain=float(ur.sum() ** 2 / (K * (ur ** 2).sum())),
                      events=events / T_SIM, bits=bits / T_SIM, minrate=float(ur.min()),
                      outage=float(np.mean(snr_db < U.GAMMA_MIN_DB)), tdec_ms=1e3 * float(np.mean(tdec)))
    return out


def agg(res):
    return {m: {k: float(np.mean([r[m][k] for r in res])) for k in res[0][m]} |
               {"rate_ci95": float(1.96 * np.std([r[m]["rate"] for r in res]) / np.sqrt(len(res)))}
            for m in METHODS}


if __name__ == "__main__":
    MC = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    Ks, speeds, tsws = [2, 4, 8], [10, 20, 30, 40], [0.1e-3, 0.2e-3, 0.5e-3, 1.0e-3]
    jobs = {}
    for K in Ks: jobs[("K", K)] = [(K, 25, 7000 + 100 * K + i, 0.2e-3) for i in range(MC)]
    for v in speeds: jobs[("v", v)] = [(4, v, 9000 + 100 * v + i, 0.2e-3) for i in range(MC)]
    for ts in tsws: jobs[("tsw", ts)] = [(4, 25, 20000 + i, ts) for i in range(MC)]
    flat = [(key, a) for key, l in jobs.items() for a in l]
    t0 = time.time()
    with Pool(os.cpu_count()) as pool:
        outs = pool.map(run_mu, [a for _, a in flat])
    print(f"{len(flat)} trials in {time.time()-t0:.0f}s")
    g = {}
    for (key, _), o in zip(flat, outs): g.setdefault(key, []).append(o)
    final = {"K": {K: agg(g[("K", K)]) for K in Ks}, "v": {v: agg(g[("v", v)]) for v in speeds},
             "tsw": {ts * 1e3: agg(g[("tsw", ts)]) for ts in tsws}, "MC": MC, "XI": XI, "T_SIM": T_SIM}
    json.dump(final, open("results/mu_results.json", "w"))
    for sec in ["K", "v", "tsw"]:
        for x, a in final[sec].items():
            print(sec, x, {m: (round(a[m]["rate"], 3), round(a[m]["events"], 3), round(a[m]["jain"], 3)) for m in METHODS})
