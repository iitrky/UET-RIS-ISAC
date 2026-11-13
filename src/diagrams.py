"""System-model diagrams (Fig. 1 scenario, Fig. 2 frame structure & processing flow)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Polygon, FancyArrowPatch
import os
os.makedirs("figures", exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 10})
INK, BLUE, RED, GREEN, GREY, ORANGE = "#222222", "#1f5fa8", "#c0392b", "#2e8b57", "#9a9a9a", "#d68910"

# ------------------------------------------------------------ Fig. 1 scenario
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.add_patch(Rectangle((-75, -14), 70, 14, fc="#d9d9d9", ec=INK, lw=1))          # buildings
ax.text(-40, -7, "building", ha="center", fontsize=8.5)
ax.add_patch(FancyBboxPatch((-4, 14.8), 14, 3.4, boxstyle="round,pad=0.2", fc="#7f8c8d", ec=INK, zorder=5))
ax.text(3, 16.5, "truck", ha="center", va="center", fontsize=7.5, color="white", zorder=6)
ax.add_patch(Rectangle((5, -14), 60, 14, fc="#d9d9d9", ec=INK, lw=1))
ax.add_patch(Rectangle((-5, -14), 10, 14, fc="#ececec", ec=INK, lw=1))
# road
ax.add_patch(Rectangle((-55, 9), 120, 10, fc="#f4f4f4", ec="none"))
for y in (9, 19):
    ax.plot([-55, 65], [y, y], color=INK, lw=1.2)
ax.plot([-55, 65], [14, 14], color=INK, lw=0.8, ls=(0, (6, 6)))
ax.text(60, 20.2, "road (2 lanes)", ha="right", fontsize=8.5)
# RIS on wall
ax.add_patch(Rectangle((-3, -0.6), 6, 1.2, fc=BLUE, ec=INK, lw=1))
for xx in (-2.2, 2.2):
    ax.add_patch(Rectangle((xx - 0.35, -0.6), 0.7, 1.2, fc=ORANGE, ec=INK, lw=0.6))
ax.annotate("Semi-passive RIS ($N$ elements)\n+ $M_s$ active sensing elements", xy=(0, -0.6), xytext=(18, -9),
            fontsize=8.5, arrowprops=dict(arrowstyle="->", lw=0.8), ha="left")
# RSU
ax.add_patch(Polygon([[-52, 38], [-48, 38], [-50, 46]], fc=INK))
ax.plot([-50, -50], [30, 38], color=INK, lw=2)
ax.text(-50, 47.5, "ISAC RSU ($N_t$ antennas)", ha="center", fontsize=8.5)
# vehicle
vx, vy = 18, 12
ax.add_patch(FancyBboxPatch((vx - 3.5, vy - 1.4), 7, 2.8, boxstyle="round,pad=0.3", fc=GREEN, ec=INK))
ax.annotate("", xy=(vx + 12, vy), xytext=(vx + 4.5, vy), arrowprops=dict(arrowstyle="-|>", lw=1.4, color=GREEN))
ax.text(vx + 8.5, vy - 3.2, "$\\mathbf{v}_k$", color=GREEN, fontsize=10, ha="center")
ax.text(vx, vy + 3.0, "vehicle $\\mathbf{p}_k$", ha="center", fontsize=8.5)
# links
ax.annotate("", xy=(-1, 0.8), xytext=(-49, 37), arrowprops=dict(arrowstyle="-|>", lw=1.6, color=BLUE))
ax.text(-27, 22, "$\\mathbf{G}$ (static LoS)", color=BLUE, rotation=-38, fontsize=9)
ax.annotate("", xy=(vx - 3, vy - 0.5), xytext=(1, 0.8), arrowprops=dict(arrowstyle="-|>", lw=1.6, color=BLUE))
ax.text(1.0, 7.2, "$\\mathbf{h}_k$", color=BLUE, fontsize=10)
ax.annotate("", xy=(2.2, 0.7), xytext=(vx - 2, vy - 1.6),
            arrowprops=dict(arrowstyle="-|>", lw=1.1, color=ORANGE, ls="--"))
ax.text(15, 2.5, "echo", color=ORANGE, fontsize=8.5)
ax.plot([-49, vx - 3], [37, vy + 1], color=RED, lw=1.2, ls=":")
ax.text(-22, 29.5, "direct path (blocked)", color=RED, fontsize=8.5, rotation=-19)
ax.plot([-48.5, -48.5, -4], [31, 3, 3], color=GREY, lw=1.1, ls="-.")
ax.text(-47, 4.2, "control link $R_c$ ($Nb$ bits per update)", color="#555555", fontsize=8)
# axes
ax.annotate("", xy=(-70, -20), xytext=(-75, -20), arrowprops=dict(arrowstyle="->", lw=0.8))
ax.text(-69, -21, "$x$", fontsize=9); ax.text(-76.5, -15.6, "$y$", fontsize=9)
ax.annotate("", xy=(-75, -15), xytext=(-75, -20), arrowprops=dict(arrowstyle="->", lw=0.8))
ax.set_xlim(-78, 67); ax.set_ylim(-22, 50); ax.set_aspect("equal"); ax.axis("off")
plt.tight_layout(); plt.savefig("figures/fig1_system_model.png", dpi=300); plt.close()

# ------------------------------------------------------------ Fig. 2 frame + flow
fig, ax = plt.subplots(figsize=(7.2, 5.2))
ax.set_xlim(0, 100); ax.set_ylim(0, 74); ax.axis("off")
# (a) frame timeline
ax.text(1, 71, "(a) Frame structure (frame length $T$)", fontsize=10, weight="bold")
x0 = 4
for i, rec in enumerate([True, False, False, True, False]):
    w = 18
    if rec:
        ax.add_patch(Rectangle((x0, 60), 3.2, 6, fc=RED, ec=INK, lw=0.8))
        ax.add_patch(Rectangle((x0 + 3.2, 60), w - 3.2, 6, fc="#cfe0f3", ec=INK, lw=0.8))
    else:
        ax.add_patch(Rectangle((x0, 60), w, 6, fc="#cfe0f3", ec=INK, lw=0.8))
    ax.text(x0 + w / 2 + (1.6 if rec else 0), 63, "ISAC data + echo", ha="center", va="center", fontsize=6.6)
    ax.text(x0 + w / 2, 57.6, f"frame $k{'+' + str(i) if i else ''}$", ha="center", fontsize=8)
    x0 += w
ax.add_patch(Rectangle((60, 51.5), 3, 2.6, fc=RED, ec=INK, lw=0.6))
ax.text(64, 52.8, "reconfiguration overhead $\\tau_c = Nb/R_c + T_{sw}$ (only when triggered)", fontsize=7.8, va="center")

# (b) flow
ax.text(1, 46, "(b) Per-frame processing at the RSU (UET-RIS)", fontsize=10, weight="bold")
def box(x, y, w, h, t, fc="#ffffff"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4", fc=fc, ec=INK, lw=0.9))
    ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=7.8)
def arr(x1, y1, x2, y2, t=None):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", lw=1))
    if t: ax.text((x1 + x2) / 2 + 0.8, (y1 + y2) / 2 + 0.8, t, fontsize=7.5)
box(1.5, 32, 18, 8, "Echo at RIS\nsensing elements: $\\mathbf{z}_k$\n(range, azimuth)", "#fdebd0")
box(24, 32, 17, 8, "EKF update\n$\\hat{\\mathbf{s}}_{k|k},\\ \\mathbf{P}_{k|k}$")
box(46, 32, 17, 8, "Predict $k{+}1$ &\nsigma points\n$\\{\\boldsymbol{\\chi}_i, w_i\\}_{i=0}^{4}$")
box(67, 32, 31, 8, "Expected rates over sigma points\n$\\bar R^{\\rm old}_{k}$ (current $\\boldsymbol{\\theta}$),\n$\\bar R^{\\rm fresh}_{k}$ (best candidate)", "#e8f1fb")
arr(19.5, 36, 23.5, 36); arr(41.5, 36, 45.5, 36); arr(63.5, 36, 67.5, 36)
box(68, 17, 29, 9, "Regret accumulation\n$A_k = A_{k-1} + [\\bar R^{\\rm fresh}_k-\\bar R^{\\rm old}_k]^+$", "#e8f1fb")
arr(82.5, 31.5, 82.5, 26.5)
ax.add_patch(Polygon([[44, 21.5], [54, 27.5], [64, 21.5], [54, 15.5]], fc="#fff7d6", ec=INK, lw=0.9))
ax.text(54, 21.5, "$A_k \\geq \\rho\\,\\bar R^{\\rm fresh}_k$\nor guard?", ha="center", va="center", fontsize=7.3)
arr(67.4, 21.5, 64.3, 21.5)
box(14, 3, 32, 9, "Lead-aware design at $\\hat{\\mathbf{p}}_{k+\\ell^\\star}$,\nsegmented candidate, $b$-bit\nquantisation, send via $R_c$", "#fadbd8")
box(60, 3, 37, 9, "Keep current $\\boldsymbol{\\theta}$ (no signalling),\nserve frame with zero overhead")
arr(50, 17.8, 38, 12.6, "yes"); arr(58, 17.8, 70, 12.6, "no")
ax.annotate("", xy=(10.5, 31.5), xytext=(13.5, 7.5),
            arrowprops=dict(arrowstyle="-|>", lw=0.9, ls="--", connectionstyle="angle,angleA=180,angleB=-90"))
ax.text(12, 20, "next frame:\nbeam sets echo\nquality (coupling)", fontsize=7.3)
plt.tight_layout(); plt.savefig("figures/fig2_frame_and_flow.png", dpi=300); plt.close()
print("ok")
