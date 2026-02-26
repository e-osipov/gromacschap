import subprocess
import os
import argparse
import shutil
import glob

import json

import numpy as np
from matplotlib import pyplot as pl

### CHAP defaults
# !!! Adjust configs for GROMACS

DOCKER_IMAGE = "CHAP:1.0"
CONTAINER_DATA = "/data"

###

# use Gemini to write some functions
# Main code comes from https://tutorials.gromacs.org/membrane-protein.html
def run_shell(command):
    """Helper to run commands inside the CHAP Docker container.
    Mounts the current working directory as /data inside the container."""
    workdir = os.path.abspath(os.getcwd())
    docker_cmd = (
        f"docker run --rm "
        f"-v {workdir}:{CONTAINER_DATA} "
        f"-w {CONTAINER_DATA} "
        f"{DOCKER_IMAGE} "
        f"{command}"
    )
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {command}\n{e}")

# 1. Provide user input
parser =  argparse.ArgumentParser()
parser.add_argument('-i',help='GHARMM-GUI output with GROMACS files. You should prepare system using CHARMM-GUI', default='input/gromacs')
# unused params. Keep just in case they needed
#parser.add_argument('-m', help='Path to MDP (params) folder')
#parser.add_argument('-f',help='Force field to use for protein')
#parser.add_argument('-w',help='Solvent force field',default='tip3p')
parser.add_argument('-o',help='Output folder name', default='run')
args = parser.parse_args()

# 2. Copy required files to run directory
os.mkdir(args.o)
shutil.copy2(f'{args.i}/step5_input.gro',f'{args.o}/')
shutil.copy2(f'{args.i}/step5_input.pdb',f'{args.o}/')
shutil.copy2(f'{args.i}/topol.top',f'{args.o}/')
shutil.copy2(f'{args.i}/index.ndx',f'{args.o}/')
shutil.copytree(f'{args.i}/toppar',f'{args.o}/toppar')
for f in glob.glob(f'{args.i}/*.mdp'):
    shutil.copy2(f,f'{args.o}/')

os.chdir(args.o)

# 3. Run minimization
run_shell("gmx grompp -f step6.0_minimization.mdp -o minimization.tpr -c step5_input.gro -r step5_input.gro -p topol.top")
run_shell("gmx mdrun -v -deffnm minimization")

# 4. (Optional) Report minimization values
#run_shell('gmx energy -f minimization.edr -o potential.xvg -xvg none')

print("Setup Complete. Ready for NVT/NPT equilibration.")

# 5. Run NVT
# 5.1 Step 1
run_shell('gmx grompp -f step6.1_equilibration.mdp -o step6.1_equilibration.tpr -c minimization.gro -r step5_input.gro -p topol.top -n index.ndx')
run_shell('gmx mdrun -v -deffnm step6.1_equilibration -nb cpu')

# 5.2 Step 2
run_shell('gmx grompp -f step6.2_equilibration.mdp -o step6.2_equilibration.tpr -c step6.1_equilibration.gro -r step5_input.gro -p topol.top -n index.ndx')
run_shell('gmx mdrun -v -deffnm step6.2_equilibration -nb cpu')

# 5.3 Eval energy
#run_shell('echo 17|gmx energy -f step6.2_equilibration_NVT_step2.edr -o  NVT_S2-temp.xvg -xvg none')

# 6. Run NPT
run_shell('gmx grompp -f step6.3_equilibration.mdp -o step6.3_equilibration.tpr -c step6.2_equilibration.gro -r step5_input.gro -p topol.top -n index.ndx')
run_shell('gmx mdrun -v -deffnm step6.3_equilibration -nb cpu')

# 6.2 Run step 2, relax constraints
run_shell('gmx grompp -f step6.4_equilibration.mdp -o step6.4_equilibration.tpr -c step6.3_equilibration.gro -r step5_input.gro -p topol.top -n index.ndx')
run_shell('gmx mdrun -v -deffnm step6.4_equilibration -nb cpu')

# 6.3 Run step 3 of NPT
run_shell('gmx grompp -f step6.5_equilibration.mdp -o step6.5_equilibration.tpr -c step6.4_equilibration.gro -r step5_input.gro -p topol.top -n index.ndx')
run_shell('gmx mdrun -v -deffnm step6.5_equilibration -nb cpu')

# 6.4 Run step 4 of NPT
run_shell('gmx grompp -f step6.6_equilibration.mdp -o step6.6_equilibration.tpr -c step6.5_equilibration.gro -r step5_input.gro -p topol.top -n index.ndx')
run_shell('gmx mdrun -v -deffnm step6.6_equilibration -nb cpu')

# 6.5 Check pressure
#run_shell('echo 18|gmx energy -f step6.6_equilibration_NPT_step4.edr -o  NPT_S4_Press.xvg -xvg none')
print("""The system's temperature should rise to the desired level (~303 K) and stay steady for the rest of the equilibration. While, pressure is the parameter which fluctuates a lot during the MD simulation, and its average should be statistically comparable with the pre-defined pressure (1 bar).""")

# 7. Production run
# Generate params
run_shell('gmx grompp -f step7_production.mdp -o step7_production.tpr -c step6.6_equilibration.gro -t step6.6_equilibration.cpt -p topol.top -n index.ndx')

run_shell('gmx mdrun -s step7_production -cpi')

# 8. Report generation
# RMSD
#run_shell('echo 4 4 |gmx rms -s step7_production.tpr -f step7_production_traj_comp_skip500.xtc -o production_rmsd.xvg -bin production.rmsd.dat -tu ns -fit rot+trans -xvg none')
# ENERGY
#run_shell('echo 13|gmx energy -f step7_production.edr -o production_Etot.xvg -xvg none')

# 9. Run CHAP
run_shell('chap -f traj_comp.xtc -s step7_production.tpr \
          -sel-pathway 1 -sel-solvent 16 -hydrophob-database zhu_2016')
# modify data and add residues if necessary

# 10. Draw graphs

with open("output.json") as data_file:
    data = json.load(data_file)


# begin new plot:
pl.figure("radius_profile")

# represent radius profile as black line:
pl.plot(                                                                        
    np.array(data["pathwayProfile"]["s"]),                                      
    np.array(data["pathwayProfile"]["radiusMean"]),                             
    "k-")  
pf = np.array(data["residueSummary"]["poreFacing"]["mean"]) > 0.5 
pl.scatter(                                                                     
    np.array(data["residueSummary"]["s"]["mean"])[pf],                          
    np.array(data["residueSummary"]["rho"]["mean"])[pf],                        
    c = np.array(data["residueSummary"]["hydrophobicity"])[pf],                 
    marker = "o",                                                               
    cmap = "BrBG_r")                                                            
pl.clim(                                                                        
    -max(abs(np.array(data["residueSummary"]["hydrophobicity"]))),              
    max(abs(np.array(data["residueSummary"]["hydrophobicity"]))))               
cbar = pl.colorbar()                                                            
cbar.ax.set_ylabel("Hydrophobicity (a.u.)")

# save plot to file:
pl.savefig("radius_profile.png")

# end figure definition:
pl.close("radius_profile")
