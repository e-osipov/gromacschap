#!/usr/bin/env python3
"""
GROMACS membrane protein MD pipeline + CHAP analysis.
Designed to run INSIDE the container (launched by run_chap.sh).

Based on: https://tutorials.gromacs.org/membrane-protein.html
"""

import subprocess
import os
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


INPUT = "/input"
OUTPUT = "/output"

# ── 1. Copy required files to output directory ───────────────────────
os.makedirs(OUTPUT, exist_ok=True)
shutil.copy2(f"{INPUT}/step5_input.gro", f"{OUTPUT}/")
shutil.copy2(f"{INPUT}/step5_input.pdb", f"{OUTPUT}/")
shutil.copy2(f"{INPUT}/topol.top", f"{OUTPUT}/")
shutil.copy2(f"{INPUT}/index.ndx", f"{OUTPUT}/")
os.makedirs(f"{OUTPUT}/toppar", exist_ok=True)
for f in glob.glob(f"{INPUT}/toppar/*"):
    shutil.copy2(f, f"{OUTPUT}/toppar/")
for f in glob.glob(f"{INPUT}/*.mdp"):
    shutil.copy2(f, f"{OUTPUT}/")

os.chdir(OUTPUT)

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

host_uid = int(os.environ.get("HOST_UID", 0))
host_gid = int(os.environ.get("HOST_GID", 0))
if host_uid:
    run_shell(f"chown -R {host_uid}:{host_gid} /output")

print("Done. Results saved in current directory.")
