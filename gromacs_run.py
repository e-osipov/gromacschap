import subprocess
import os
import argparse
import shutil

# use Gemini to write some functions
# Main code comes from https://tutorials.gromacs.org/membrane-protein.html
def run_gmx(command):
    """Helper to run GROMACS commands via shell"""
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {command}\n{e}")

# 1. Provide user input
parser =  argparse.ArgumentParser()
parser.add_argument('-i',help='GHARMM-GUI output with GROMACS files. You should prepare system using CHARMM-GUI', default='input/gromacs')
parser.add_argument('-f',help='Force field to use for protein')
parser.add_argument('-w',help='Solvent force field',default='tip3p')
parser.add_argument('-o',help='Output folder name', default='run')
args = parser.parse_args()

# 2. Copy required files to run directory
shutil.copy2(f'{args.i}/step5_input.gro',args.o)
shutil.copy2(f'{args.i}/step5_input.pdb',args.o)
shutil.copy2(f'{args.i}/topol.top',args.o)
shutil.copy2(f'{args.i}index.ndx',args.o)
shutil.copy2(f'{args.i}/toppar',args.o)
shutil.copy2('tutorial/mdp/*.mdp',args.o)

os.chdir(args.o)

# 3. Run minimization
run_gmx("gmx grompp -f step6.0_minimization.mdp -o minimization.tpr -c step5_input.gro -r step5_input.gro -p topol.top")
run_gmx("gmx mdrun -v -deffnm minimization")

# 4. Report minimization values
run_gmx('gmx energy -f minimization.edr -o potential.xvg -xvg none')

print("Setup Complete. Ready for NVT/NPT equilibration.")

# 5. Run NVT
# 5.1 Step 1
run_gmx('gmx grompp -f step6.1_equilibration_NVT_step1.mdp -o step6.1_equilibration_NVT_step1.tpr -c minimization.gro -r step5_input.gro -p topol.top -n index.ndx')
run_gmx('gmx mdrun -v -deffnm step6.1_equilibration_NVT_step1')

# 5.2 Step 2
run_gmx('gmx grompp -f step6.2_equilibration_NVT_step2.mdp -o step6.2_equilibration_NVT_step2.tpr -c step6.1_equilibration_NVT_step1.gro -r step5_input.gro -p topol.top -n index.ndx')
run_gmx('gmx mdrun -v -deffnm step6.2_equilibration_NVT_step2')

# 5.3 Eval energy
run_gmx('gmx energy -f step6.2_equilibration_NVT_step2.edr -o  NVT_S2-temp.xvg -xvg none')

# 6. Run NPT
run_gmx('gmx grompp -f step6.3_equilibration_NPT_step1.mdp -o step6.3_equilibration_NPT_step1.tpr -c step6.2_equilibration_NVT_step2.gro -r step5_input.gro -p topol.top -n index.ndx')
run_gmx('gmx mdrun -v -deffnm step6.3_equilibration_NPT_step1')

# 6.2 Run step 2, reduce constraints
run_gmx('gmx grompp -f step6.4_equilibration_NPT_step2.mdp -o step6.4_equilibration_NPT_step2.tpr -c step6.3_equilibration_NPT_step1.gro -r step5_input.gro -p topol.top -n index.ndx')
run_gmx('gmx mdrun -v -deffnm step6.4_equilibration_NPT_step2')

# 6.3 Run step 3 of NPT
run_gmx('gmx grompp -f step6.5_equilibration_NPT_step3.mdp -o step6.5_equilibration_NPT_step3.tpr -c step6.4_equilibration_NPT_step2.gro -r step5_input.gro -p topol.top -n index.ndx')
run_gmx('gmx mdrun -v -deffnm step6.5_equilibration_NPT_step3')

# 6.4 Run step 4 of NPT
run_gmx('gmx grompp -f step6.6_equilibration_NPT_step4.mdp -o step6.6_equilibration_NPT_step4.tpr -c step6.5_equilibration_NPT_step3.gro -r step5_input.gro -p topol.top -n index.ndx')
run_gmx('gmx mdrun -v -deffnm step6.6_equilibration_NPT_step4')

# 6.5 Check pressure
run_gmx('echo 18|gmx energy -f step6.6_equilibration_NPT_step4.edr -o  NPT_S4_Press.xvg -xvg none')
print("""The system's temperature should rise to the desired level (~303 K) and stay steady for the rest of the equilibration. While, pressure is the parameter which fluctuates a lot during the MD simulation, and its average should be statistically comparable with the pre-defined pressure (1 bar).""")

# 7. Production run
# Generate params
run_gmx('gmx grompp -f step7_production_revised.mdp -o step7_production.tpr -c step6.6_equilibration_NPT_step4.gro -t step6.6_equilibration_NPT_step4.cpt -p topol.top -n index.ndx')

run_gmx('gmx mdrun -s step7_production -cpi')

# 8. Report generation
# RMSD
run_gmx('echo 4 4 |gmx rms -s step7_production.tpr -f step7_production_traj_comp_skip500.xtc -o production_rmsd.xvg -bin production.rmsd.dat -tu ns -fit rot+trans -xvg none')
# ENERGY
run_gmx('echo 13|gmx energy -f step7_production.edr -o production_Etot.xvg -xvg none')
