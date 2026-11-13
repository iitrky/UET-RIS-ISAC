# UET-RIS: Regret Triggered RIS Reconfiguration for ISAC Vehicular Networks

Official simulation code for the paper

> G. Arul Dalton and S. V. J. Mani, "Regret Triggered RIS Reconfiguration for ISAC Vehicular Networks," manuscript submitted for publication, 2026.

The repository reproduces **every figure and every number** reported in the paper: Figs. 1–10 and the numerical results in Tables 5–6.

---

## 1. What the code does

A reconfigurable intelligent surface (RIS) relays a millimetre-wave (28 GHz) link from an ISAC roadside unit (RSU) to vehicles whose direct path is blocked. Every RIS update must be sent over a rate-limited control link and then applied, which costs part of the frame:

```
tau_c = N*b / R_c + T_sw
```

The code answers the question **when should the RIS be reconfigured?** It implements:

| Component | Paper section | File |
|---|---|---|
| Exact spherical-wavefront RSU–RIS–vehicle channel, *b*-bit phases | Sec. 3.2 | `src/uet_ris_sim.py` (`RIS`) |
| Echo-based EKF tracking with sensing–illumination coupling | Sec. 3.4 | `src/uet_ris_sim.py` (`h_meas`, `H_jac`, `sens_std`) |
| **UET-RIS**: sigma-point expected rate, regret (ski-rental) trigger, lead-compensated focusing | Sec. 5 | `src/uet_ris_sim.py` (`run_trial`, method `uet`) |
| Hold-time scaling law validation, *H\** ∝ *v*^(−2/3) | Sec. 5.5, 8.6 | `src/scaling.py` |
| **Multi-vehicle batched trigger** (online joint replenishment) | Sec. 6 | `src/mu_sim.py` (method `jrp`) |
| Baselines (periodic M = 1/5/10/20, mean-only ablation, independent triggers, genie) | Sec. 8.1, 8.7 | both simulators |

Implementation details, including the mapping from equations to code, are in **[`docs/IMPLEMENTATION.md`](docs/IMPLEMENTATION.md)**.

---

## 2. Repository structure

```
UET-RIS-ISAC/
├── README.md               this file
├── LICENSE                 MIT licence
├── CITATION.cff            citation metadata (GitHub "Cite this repository")
├── requirements.txt        Python dependencies
├── reproduce.py            one-command reproduction (full or --quick)
├── run_all.sh              same as reproduce.py, for Linux/macOS shells
├── src/
│   ├── uet_ris_sim.py      single-vehicle simulator (Figs. 3–8)
│   ├── mu_sim.py           multi-vehicle simulator (Figs. 9–10, Table 6)
│   ├── plots.py            plots Figs. 3–7 from results/results.json
│   ├── scaling.py          plots Fig. 8 (hold-time law)
│   ├── mu_plots.py         plots Figs. 9–10 from results/mu_results.json
│   └── diagrams.py         draws the system-model diagrams (Figs. 1–2)
├── results/
│   ├── results.json        saved single-vehicle results used in the paper
│   └── mu_results.json     saved multi-vehicle results used in the paper
├── figures/                all paper figures (PNG, 300 dpi)
└── docs/
    └── IMPLEMENTATION.md   detailed implementation notes
```

---

## 3. Installation

Requires **Python 3.9 or later**. Only NumPy and Matplotlib are needed.

```bash
git clone https://github.com/iitrky/UET-RIS-ISAC.git
cd UET-RIS-ISAC
python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows:
# .venv\Scripts\activate
pip install -r requirements.txt
```

---

## 4. Reproducing the results

### 4.1 Quick check (about 1–2 minutes)

Runs both simulators with 3 Monte Carlo trials per point in a separate folder (`quick_check/`). It does **not** overwrite the published results.

```bash
python reproduce.py --quick
```

### 4.2 Full reproduction (about 15–20 minutes on a 2-core CPU)

```bash
python reproduce.py
```

This runs, in order:

| Step | Command | Output | Paper item |
|---|---|---|---|
| 1 | `python src/diagrams.py` | `figures/fig1_*.png`, `fig2_*.png` | Figs. 1–2 |
| 2 | `python src/uet_ris_sim.py 40` | `results/results.json` | data for Figs. 3–8, Table 5 |
| 3 | `python src/plots.py` | `figures/fig3`–`fig7` | Figs. 3–7 |
| 4 | `python src/scaling.py` | `figures/fig8_scaling_law.png` | Fig. 8 |
| 5 | `python src/mu_sim.py 30` | `results/mu_results.json` | data for Figs. 9–10, Table 6 |
| 6 | `python src/mu_plots.py` | `figures/fig9`, `fig10` | Figs. 9–10 |

You can also run any step on its own from the repository root, for example `python src/uet_ris_sim.py 40`. The number after the script name is the number of Monte Carlo trials per point.

All random seeds are fixed, so a full run reproduces the saved JSON files. Floating-point results may differ in the last digits across NumPy versions and CPUs.

---

## 5. Main parameters (Table 4 of the paper)

| Parameter | Value | Variable in `src/uet_ris_sim.py` |
|---|---|---|
| Carrier / bandwidth | 28 GHz / 50 MHz | `FC`, `BW` |
| Noise figure | 7 dB (noise −90 dBm) | `NF_DB` |
| RSU power / antennas | 33 dBm / 32 | `PT_W`, `NT` |
| RIS size | 32×32 = 1024 (default) | argument `nx` |
| Phase resolution | 2 bits | `BITS` |
| Frame length | 10 ms | `T_FRAME` |
| Control-link rate / switching time | 1 Mb/s / 0.2 ms | `R_CTRL`, `T_SW` |
| Angle / range accuracy at 20 m | 1° / 0.1 m | `SIG_TH0`, `SIG_R0`, `D0` |
| EKF process noise | 1 m/s² | `SIG_A` |
| Reliability threshold, max. lead | 3 dB, 5 frames | `GAMMA_MIN_DB`, `LEAD_MAX` |
| Batching threshold ξ (multi-vehicle) | 0.5 | `XI` in `src/mu_sim.py` |

---

## 6. Method names used in the code

| Code name | Paper label |
|---|---|
| `genie` | Genie bound (true position, zero overhead) |
| `periodic1`, `periodic5`, `periodic10`, `periodic20` | Periodic, M = 1 / 5 / 10 / 20 |
| `et_mean` | Proposed trigger, mean-only (ablation) |
| `uet` | **Proposed UET-RIS** |
| `random` | Random RIS phases (reference only) |
| `indep` (multi-vehicle) | Independent regret triggers |
| `jrp` (multi-vehicle) | **Proposed JRP-batched UET-RIS** |

---

## 7. Citation

If you use this code, please cite the paper (the BibTeX entry will be updated once the paper is published):

```bibtex
@article{dalton2026uetris,
  author  = {Arul Dalton, G. and Mani, S. V. J.},
  title   = {Regret Triggered {RIS} Reconfiguration for {ISAC} Vehicular Networks},
  journal = {Manuscript submitted for publication},
  year    = {2026}
}
```

---

## 8. Licence and contact

Released under the MIT licence (see `LICENSE`).

Corresponding author: **G. Arul Dalton**, Department of Computer Science and Engineering, Saveetha Engineering College, Chennai, India. aruldalton@saveetha.ac.in · ORCID [0000-0002-5981-7176](https://orcid.org/0000-0002-5981-7176)

Co-author: **S. V. J. Mani**, Department of Computer Science Engineering and Cyber Security, Saveetha Engineering College, Chennai, India. manisvj@saveetha.ac.in · ORCID [0009-0008-0957-6301](https://orcid.org/0009-0008-0957-6301)
