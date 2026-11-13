"""Figs. 9-10: multi-vehicle results from results/mu_results.json"""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 10, "axes.grid": True, "grid.alpha": 0.3})
R = json.load(open("results/mu_results.json"))
STY = {"genie": ("Genie focusing (true positions, zero overhead)", "k", "--", ""),
       "periodic1": ("Periodic, M=1", "#7f7f7f", "-", "v"), "periodic5": ("Periodic, M=5", "#1f77b4", "-", "s"),
       "indep": ("Independent regret triggers", "#ff7f0e", "-.", "x"), "jrp": ("Proposed JRP-batched UET-RIS", "#d62728", "-", "o")}
def two(sec, xs, xl, fname, xt=None):
    fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.9))
    for m, (lab, c, ls, mk) in STY.items():
        X = xt if xt else xs
        ax[0].plot(X, [R[sec][str(x)][m]["rate"] for x in xs], color=c, ls=ls, marker=mk, label=lab)
        if m != "genie":
            ax[1].plot(X, [R[sec][str(x)][m]["events"] for x in xs], color=c, ls=ls, marker=mk, label=lab)
    ax[0].set_ylabel("Avg. per-vehicle effective SE (bit/s/Hz)"); ax[1].set_ylabel("Switching events per frame")
    for a in ax: a.set_xlabel(xl)
    ax[0].set_title("(a)", fontsize=10); ax[1].set_title("(b)", fontsize=10)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=7.5, frameon=False)
    plt.tight_layout(rect=(0, 0.1, 1, 1)); plt.savefig(f"figures/{fname}", dpi=300); plt.close()
two("K", [2, 4, 8], "Number of vehicles K (sub-arrays)", "fig9_mu_vs_K.png")
two("tsw", [0.1, 0.2, 0.5, 1.0], "RIS switching/settling time T_sw (ms)", "fig10_mu_vs_Tsw.png")
print("ok")
