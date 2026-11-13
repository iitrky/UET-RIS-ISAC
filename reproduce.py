"""
One-command reproduction of all figures and results in the paper.

    python reproduce.py            # full run (about 15-20 min on 2 cores), overwrites results/ and figures/
    python reproduce.py --quick    # smoke test with 3 trials per point, written to quick_check/ only
"""
import os, sys, shutil, subprocess, time

ROOT = os.path.dirname(os.path.abspath(__file__))
STEPS = [("diagrams.py", []), ("uet_ris_sim.py", ["{MC1}"]), ("plots.py", []),
         ("scaling.py", []), ("mu_sim.py", ["{MC2}"]), ("mu_plots.py", [])]


def main():
    quick = "--quick" in sys.argv
    if quick:
        work = os.path.join(ROOT, "quick_check")
        shutil.rmtree(work, ignore_errors=True)
        shutil.copytree(os.path.join(ROOT, "src"), os.path.join(work, "src"))
        mc1, mc2 = "3", "3"
    else:
        work, mc1, mc2 = ROOT, "40", "30"
    os.makedirs(os.path.join(work, "results"), exist_ok=True)
    os.makedirs(os.path.join(work, "figures"), exist_ok=True)
    t_all = time.time()
    for script, args in STEPS:
        args = [a.replace("{MC1}", mc1).replace("{MC2}", mc2) for a in args]
        cmd = [sys.executable, os.path.join("src", script)] + args
        print(f"\n=== {' '.join(cmd)}  (in {work})", flush=True)
        t0 = time.time()
        subprocess.run(cmd, cwd=work, check=True)
        print(f"    done in {time.time() - t0:.0f} s", flush=True)
    print(f"\nAll steps finished in {time.time() - t_all:.0f} s. Figures: {os.path.join(work, 'figures')}")


if __name__ == "__main__":
    main()
