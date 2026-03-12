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


def step_done(*output_files):
    """Return True if all output files for a step already exist."""
    return all(os.path.exists(f) for f in output_files)


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
if step_done(f"{OUTPUT}/minimization.gro"):
    print("Minimization results found, skipping.")
else:
    run_shell(
        "gmx grompp -f step6.0_minimization.mdp -o minimization.tpr "
        "-c step5_input.gro -r step5_input.gro -p topol.top"
    )
    run_shell("gmx mdrun -v -deffnm minimization")

print("Minimization complete. Starting equilibration...")

# ── 4. NVT equilibration (steps 1–2) ────────────────────────────────
if step_done(f"{OUTPUT}/step6.1_equilibration.gro"):
    print("step6.1 equilibration results found, skipping.")
else:
    run_shell(
        "gmx grompp -f step6.1_equilibration.mdp -o step6.1_equilibration.tpr "
        "-c minimization.gro -r step5_input.gro -p topol.top -n index.ndx"
    )
    run_shell("gmx mdrun -v -deffnm step6.1_equilibration ")

if step_done(f"{OUTPUT}/step6.2_equilibration.gro"):
    print("step6.2 equilibration results found, skipping.")
else:
    run_shell(
        "gmx grompp -f step6.2_equilibration.mdp -o step6.2_equilibration.tpr "
        "-c step6.1_equilibration.gro -r step5_input.gro -p topol.top -n index.ndx"
    )
    run_shell("gmx mdrun -v -deffnm step6.2_equilibration ")

# ── 5. NPT equilibration (steps 3–6) ────────────────────────────────
for step in range(3, 7):
    prev = step - 1 if step > 3 else 2
    prev_name = f"step6.{prev}_equilibration"
    curr_name = f"step6.{step}_equilibration"
    if step_done(f"{OUTPUT}/{curr_name}.gro"):
        print(f"{curr_name} results found, skipping.")
        continue
    run_shell(
        f"gmx grompp -f {curr_name}.mdp -o {curr_name}.tpr "
        f"-c {prev_name}.gro -r step5_input.gro -p topol.top -n index.ndx"
    )
    run_shell(f"gmx mdrun -v -deffnm {curr_name} ")

print(
    "Equilibration complete. Temperature should be ~303 K and "
    "average pressure ~1 bar (large fluctuations are normal)."
)

# ── 6. Production run ───────────────────────────────────────────────
if step_done(f"{OUTPUT}/traj_comp.xtc"):
    print("Production run results found, skipping.")
else:
    run_shell(
        "gmx grompp -f step7_production.mdp -o step7_production.tpr "
        "-c step6.6_equilibration.gro -t step6.6_equilibration.cpt "
        "-p topol.top -n index.ndx"
    )
    run_shell("gmx mdrun -s step7_production -cpi")

# Convert output to multiframe PDB
run_shell("gmx trjconv -s step7_production.tpr -f traj_comp.xtc -o whole.xtc -pbc whole")
run_shell("gmx trjconv -s step7_production.tpr -f whole.xtc -o clean.xtc -center -pbc mol -ur compact")
run_shell("echo 0 | gmx trjconv -s step7_production.tpr -f clean.xtc -o step7_production.pdb -dt 100")

# ── 7. CHAP analysis ────────────────────────────────────────────────
if step_done(f"{OUTPUT}/output.json"):
    print("CHAP analysis results found, skipping.")
else:
    run_shell(
        "chap -f traj_comp.xtc -s step7_production.tpr "
        "-sel-pathway 1 -sel-solvent 16 -hydrophob-database zhu_2016 "
        "-hydrophob-fallback 0.0"
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
hydro_lim = max(abs(hydro))

# Radius profile with variance bands
pl.figure("radius_profile_full")
pl.plot(
    np.array(data["pathwayProfile"]["s"]),
    np.array(data["pathwayProfile"]["radiusMean"]),
    "k-",
)
pl.fill_between(
    np.array(data["pathwayProfile"]["s"]),
    np.array(data["pathwayProfile"]["radiusMin"]),
    np.array(data["pathwayProfile"]["radiusMax"]),
    facecolor="#000000",
    alpha=0.1,
)
radius_sd = np.array(data["pathwayProfile"]["radiusSd"])
pl.fill_between(
    np.array(data["pathwayProfile"]["s"]),
    np.array(data["pathwayProfile"]["radiusMean"]) - radius_sd,
    np.array(data["pathwayProfile"]["radiusMean"]) + radius_sd,
    facecolor="#000000",
    alpha=0.2,
)
pl.scatter(
    np.array(data["residueSummary"]["s"]["mean"])[pf],
    np.array(data["residueSummary"]["rho"]["mean"])[pf],
    c=hydro[pf],
    marker="o",
    cmap="BrBG_r",
)
pl.clim(-hydro_lim, hydro_lim)
cbar = pl.colorbar()
cbar.ax.set_ylabel("Hydrophobicity (a.u.)")
pl.margins(x=0)
pl.title("Time-Averaged Radius Profile")
pl.xlabel("s (nm)")
pl.ylabel("R (nm)")
pl.savefig("time_averaged_radius_profile.png", dpi=300)
pl.close("radius_profile_full")

# Hydrophobicity profile
pl.figure("hydrophobicity_profile")
pl.plot(
    np.array(data["pathwayProfile"]["s"]),
    np.array(data["pathwayProfile"]["pfHydrophobicityMean"]),
    "k-",
)
pl.fill_between(
    np.array(data["pathwayProfile"]["s"]),
    np.array(data["pathwayProfile"]["pfHydrophobicityMin"]),
    np.array(data["pathwayProfile"]["pfHydrophobicityMax"]),
    facecolor="#000000",
    alpha=0.1,
)
hydrophobicity_sd = np.array(data["pathwayProfile"]["pfHydrophobicitySd"])
pl.fill_between(
    np.array(data["pathwayProfile"]["s"]),
    np.array(data["pathwayProfile"]["pfHydrophobicityMean"]) - hydrophobicity_sd,
    np.array(data["pathwayProfile"]["pfHydrophobicityMean"]) + hydrophobicity_sd,
    facecolor="#000000",
    alpha=0.2,
)
pl.scatter(
    np.array(data["residueSummary"]["s"]["mean"])[pf],
    hydro[pf],
    c=hydro[pf],
    marker="o",
    cmap="BrBG_r",
)
pl.clim(-hydro_lim, hydro_lim)
cbar = pl.colorbar()
cbar.ax.set_ylabel("Hydrophobicity (a.u.)")
pl.margins(x=0)
pl.title("Time-Averaged Hydrophobicity Profile")
pl.xlabel("s (nm)")
pl.ylabel("H (a.u.)")
pl.savefig("time_averaged_hydrophobicity_profile.png", dpi=300)
pl.close("hydrophobicity_profile")

# Solvent number density profile
pl.figure("density_profile")
pl.axhline(y=33.3679, linestyle="dashed")
pl.plot(
    np.array(data["pathwayProfile"]["s"]),
    np.array(data["pathwayProfile"]["densityMean"]),
    "k-",
)
pl.fill_between(
    np.array(data["pathwayProfile"]["s"]),
    np.array(data["pathwayProfile"]["densityMin"]),
    np.array(data["pathwayProfile"]["densityMax"]),
    facecolor="#000000",
    alpha=0.1,
)
density_sd = np.array(data["pathwayProfile"]["densitySd"])
pl.fill_between(
    np.array(data["pathwayProfile"]["s"]),
    np.array(data["pathwayProfile"]["densityMean"]) - density_sd,
    np.array(data["pathwayProfile"]["densityMean"]) + density_sd,
    facecolor="#000000",
    alpha=0.2,
)
pl.margins(x=0)
pl.title("Time-Averaged Solvent Number Density Profile")
pl.xlabel("s (nm)")
pl.ylabel("n (nm$^{-3}$)")
pl.savefig("time_averaged_solvent_number_density_profile.png", dpi=300)
pl.close("density_profile")

# Free energy profile
pl.figure("energy_profile")
pl.plot(
    np.array(data["pathwayProfile"]["s"]),
    np.array(data["pathwayProfile"]["energyMean"]),
    "k-",
)
energy_sd = np.array(data["pathwayProfile"]["energySd"])
pl.fill_between(
    np.array(data["pathwayProfile"]["s"]),
    np.array(data["pathwayProfile"]["energyMean"]) - energy_sd,
    np.array(data["pathwayProfile"]["energyMean"]) + energy_sd,
    facecolor="#000000",
    alpha=0.2,
)
pl.margins(x=0)
pl.title("Time-Averaged Free Energy Profile")
pl.xlabel("s (nm)")
pl.ylabel("G (kT)")
pl.savefig("time_averaged_free_energy_profile.png", dpi=300)
pl.close("energy_profile")

# ── 10. Save plotted data as CSV ────────────────────────────────────
pp = data["pathwayProfile"]
s = np.array(pp["s"])
pathway_data = np.column_stack([
    s,
    np.array(pp["radiusMean"]),
    np.array(pp["radiusMin"]),
    np.array(pp["radiusMax"]),
    np.array(pp["radiusSd"]),
    np.array(pp["pfHydrophobicityMean"]),
    np.array(pp["pfHydrophobicityMin"]),
    np.array(pp["pfHydrophobicityMax"]),
    np.array(pp["pfHydrophobicitySd"]),
    np.array(pp["densityMean"]),
    np.array(pp["densityMin"]),
    np.array(pp["densityMax"]),
    np.array(pp["densitySd"]),
    np.array(pp["energyMean"]),
    np.array(pp["energySd"]),
])
np.savetxt(
    "pathway_profile.csv",
    pathway_data,
    delimiter=",",
    header=(
        "s_nm,radiusMean_nm,radiusMin_nm,radiusMax_nm,radiusSd_nm,"
        "pfHydrophobicityMean,pfHydrophobicityMin,pfHydrophobicityMax,pfHydrophobicitySd,"
        "densityMean_nm-3,densityMin_nm-3,densityMax_nm-3,densitySd_nm-3,"
        "energyMean_kT,energySd_kT"
    ),
    comments="",
)

rs = data["residueSummary"]
resi = np.array(rs['id'])[pf]
res = np.array(rs['name'])[pf]
residue_labels = np.array([f"{name}{rid}" for name, rid in zip(res, resi)])

residue_data = np.column_stack([
    np.array(rs["s"]["mean"])[pf],
    np.array(rs["rho"]["mean"])[pf],
    hydro[pf],
    residue_labels
])
np.savetxt(
    "residue_summary_pore_facing.csv",
    residue_data,
    delimiter=",",
    header="s_mean_nm,rho_mean_nm,hydrophobicity,residue labels",
    comments="",
)

print("Done. Results saved in current directory.")
