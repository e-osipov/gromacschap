#!/usr/bin/env python3
"""
GROMACS membrane protein MD pipeline + CHAP analysis.
Designed to run INSIDE the container (launched by run_chap.sh).

Based on: https://tutorials.gromacs.org/membrane-protein.html
"""

import subprocess
import os
import argparse
import shutil
import glob
import json

import numpy as np
from matplotlib import pyplot as pl


def run_shell(command):
    """Run a shell command directly. Exits on failure."""
    print(f">>> {command}")
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {command}\n{e}")
        raise SystemExit(1)


# ── 1. Parse user input ─────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Run GROMACS equilibration + production + CHAP analysis"
)
parser.add_argument(
    "-i",
    help="Input directory with CHARMM-GUI/GROMACS files (mounted as /input in container)",
    default="/input",
)
parser.add_argument("-o", help="Output directory (mounted as /output in container)", default="/output")
args = parser.parse_args()

# ── 2. Copy required files to run directory ──────────────────────────
os.makedirs(args.o, exist_ok=True)
shutil.copy2(f"{args.i}/step5_input.gro", f"{args.o}/")
shutil.copy2(f"{args.i}/step5_input.pdb", f"{args.o}/")
shutil.copy2(f"{args.i}/topol.top", f"{args.o}/")
shutil.copy2(f"{args.i}/index.ndx", f"{args.o}/")
shutil.copytree(f"{args.i}/toppar", f"{args.o}/toppar", dirs_exist_ok=True)
for f in glob.glob(f"{args.i}/*.mdp"):
    shutil.copy2(f, f"{args.o}/")

os.chdir(args.o)

# ── 3. Minimization ─────────────────────────────────────────────────
run_shell(
    "gmx grompp -f step6.0_minimization.mdp -o minimization.tpr "
    "-c step5_input.gro -r step5_input.gro -p topol.top"
)
run_shell("gmx mdrun -v -deffnm minimization")

print("Minimization complete. Starting equilibration...")

# ── 4. NVT equilibration (steps 1–2) ────────────────────────────────
run_shell(
    "gmx grompp -f step6.1_equilibration.mdp -o step6.1_equilibration.tpr "
    "-c minimization.gro -r step5_input.gro -p topol.top -n index.ndx"
)
run_shell("gmx mdrun -v -deffnm step6.1_equilibration -nb cpu")

run_shell(
    "gmx grompp -f step6.2_equilibration.mdp -o step6.2_equilibration.tpr "
    "-c step6.1_equilibration.gro -r step5_input.gro -p topol.top -n index.ndx"
)
run_shell("gmx mdrun -v -deffnm step6.2_equilibration -nb cpu")

# ── 5. NPT equilibration (steps 3–6) ────────────────────────────────
for step in range(3, 7):
    prev = step - 1 if step > 3 else 2
    prev_name = f"step6.{prev}_equilibration"
    curr_name = f"step6.{step}_equilibration"
    run_shell(
        f"gmx grompp -f {curr_name}.mdp -o {curr_name}.tpr "
        f"-c {prev_name}.gro -r step5_input.gro -p topol.top -n index.ndx"
    )
    run_shell(f"gmx mdrun -v -deffnm {curr_name} -nb cpu")

print(
    "Equilibration complete. Temperature should be ~303 K and "
    "average pressure ~1 bar (large fluctuations are normal)."
)

# ── 6. Production run ───────────────────────────────────────────────
run_shell(
    "gmx grompp -f step7_production.mdp -o step7_production.tpr "
    "-c step6.6_equilibration.gro -t step6.6_equilibration.cpt "
    "-p topol.top -n index.ndx"
)
run_shell("gmx mdrun -s step7_production -cpi")

# ── 7. CHAP analysis ────────────────────────────────────────────────
run_shell(
    "chap -f traj_comp.xtc -s step7_production.tpr "
    "-sel-pathway 1 -sel-solvent 16 -hydrophob-database zhu_2016"
)

# ── 8. Plot radius profile ──────────────────────────────────────────
with open("output.json") as f:
    data = json.load(f)

pl.figure("radius_profile")
pl.plot(
    np.array(data["pathwayProfile"]["s"]),
    np.array(data["pathwayProfile"]["radiusMean"]),
    "k-",
)

pf = np.array(data["residueSummary"]["poreFacing"]["mean"]) > 0.5
hydro = np.array(data["residueSummary"]["hydrophobicity"])

pl.scatter(
    np.array(data["residueSummary"]["s"]["mean"])[pf],
    np.array(data["residueSummary"]["rho"]["mean"])[pf],
    c=hydro[pf],
    marker="o",
    cmap="BrBG_r",
)
pl.clim(-max(abs(hydro)), max(abs(hydro)))
cbar = pl.colorbar()
cbar.ax.set_ylabel("Hydrophobicity (a.u.)")
pl.xlabel("s (nm)")
pl.ylabel("Radius (nm)")
pl.savefig("radius_profile.png", dpi=150)
pl.close("radius_profile")

print("Done. Results saved in current directory.")
