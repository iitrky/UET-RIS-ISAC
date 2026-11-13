"""Result figures from results/results.json"""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 10, "axes.grid": True,
                     "grid.alpha": 0.3, "legend.fontsize": 8})
R = json.load(open("results/results.json"))
STY = {"genie": ("Genie bound (perfect position, zero overhead)", "k", "--", ""),
       "periodic1": ("Periodic, M=1 [predictive BF]", "#7f7f7f", "-", "v"),
       "periodic5": ("Periodic, M=5", "#1f77b4", "-", "s"),
       "periodic10": ("Periodic, M=10", "#2ca02c", "-", "^"),
       "periodic20": ("Periodic, M=20", "#9467bd", "-", "D"),
       "et_mean": ("Proposed trigger, mean-only (ablation)", "#ff7f0e", "-.", "x"),
       "uet": ("Proposed UET-RIS", "#d62728", "-", "o")}
ORDER = ["genie", "periodic1", "periodic5", "periodic10", "periodic20", "et_mean", "uet"]


def sweep(key, xs, field, xlabel, ylabel, fname, xt=None, logy=False, methods=ORDER, legend_loc="best"):
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    for m in methods:
        lab, c, ls, mk = STY[m]
        y = [R[key][str(x)][m][field] for x in xs]
        ax.plot(xt if xt else xs, y, color=c, ls=ls, marker=mk, ms=5, lw=1.6 if m == "uet" else 1.2, label=lab)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if logy: ax.set_yscale("log")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=7, frameon=False)
    plt.tight_layout(); plt.savefig(f"figures/{fname}", dpi=300); plt.close()


speeds = [10, 15, 20, 25, 30, 35, 40]
sizes = [8, 12, 16, 24, 32]
sens = [0.5, 1.0, 2.0, 3.0, 4.0]
sweep("speed", speeds, "rate", "Vehicle speed (m/s)", "Average effective rate (bit/s/Hz)", "fig3_rate_vs_speed.png")
sweep("size", sizes, "rate", "Number of RIS elements N", "Average effective rate (bit/s/Hz)",
      "fig4_rate_vs_N.png", xt=[s * s for s in sizes])
sweep("speed", speeds, "recfg", "Vehicle speed (m/s)", "Reconfiguration frequency (per frame)",
      "fig5_recfg_vs_speed.png", methods=["periodic1", "periodic5", "periodic10", "periodic20", "et_mean", "uet"], logy=True)
# sensing accuracy: rate and outage
fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.6))
for m in ["genie", "periodic5", "periodic10", "et_mean", "uet"]:
    lab, c, ls, mk = STY[m]
    axs[0].plot(sens, [R["sens"][str(t)][m]["rate"] for t in sens], color=c, ls=ls, marker=mk, label=lab)
    if m != "genie":
        axs[1].plot(sens, [max(R["sens"][str(t)][m]["outage"], 1e-4) for t in sens], color=c, ls=ls, marker=mk, label=lab)
axs[0].set_xlabel(r"Angle accuracy at 20 m, $\sigma_{\theta,0}$ (deg)"); axs[0].set_ylabel("Average effective rate (bit/s/Hz)")
axs[1].set_xlabel(r"Angle accuracy at 20 m, $\sigma_{\theta,0}$ (deg)"); axs[1].set_ylabel("Outage probability (SNR < 3 dB)")
axs[1].set_yscale("log"); axs[0].legend(fontsize=7)
axs[0].set_title("(a)", fontsize=10); axs[1].set_title("(b)", fontsize=10)
plt.tight_layout(); plt.savefig("figures/fig6_sensing_accuracy.png", dpi=300); plt.close()
# trace
tr = R["trace"]; traj = np.array(tr["traj"]); x = traj[:, 0]
fig, axs = plt.subplots(2, 1, figsize=(6.4, 5.0), sharex=True)
for m in ["genie", "periodic1", "periodic5", "uet"]:
    lab, c, ls, mk = STY[m]
    axs[0].plot(x, tr[m]["trace_snr"], color=c, ls=ls, lw=1.0, label=lab)
snr_u = np.array(tr["uet"]["trace_snr"]); rr = np.array(tr["uet"]["trace_rate"])
# reconfiguration instants = frames where effective rate < log2(1+snr) (overhead charged)
idx = np.where(rr < np.log2(1 + 10 ** (snr_u / 10)) - 1e-9)[0]
axs[0].plot(x[idx], snr_u[idx], "o", color=STY["uet"][1], ms=3.5, mfc="none", label="UET-RIS reconfiguration instants")
axs[0].set_ylabel("Received SNR (dB)"); axs[0].legend(fontsize=7, loc="lower center", ncol=2)
axs[1].plot(x, tr["uet"]["trace_err"], color=STY["uet"][1], lw=1, label="UET-RIS")
axs[1].plot(x, tr["periodic5"]["trace_err"], color=STY["periodic5"][1], lw=1, label="Periodic, M=5")
axs[1].set_ylabel("Position error (m)"); axs[1].set_xlabel("Vehicle position along road x (m)"); axs[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig("figures/fig7_trace.png", dpi=300); plt.close()
print("plots ok")
