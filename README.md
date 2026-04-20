# gromacschap
Gromacs and CHAP pipeline  automated pipeline

# Installation
- Install Podman >4.2 with GPU support (ask Claude for specific commands for your distro)
- Fetch code using `git clone`
- Adjust `Dockerfile` for your architecture:
    - first, find out what is your GPU architecture. Code for Nvidia GPU only: `nvidia-smi --query-gpu=gpu_name,compute_cap --format=csv`. For example, RTX4000 has compute 8.9
    - Multipy number from previous step by 10 and change "89" in `Dockerfile`: `-DGMX_CUDA_TARGET_COMPUTE="89" ` to your value
- Build image: `podman build  --format docker -t chap:1.0 -t chap:latest .`
- (Optional) Add location of the repo to `$PATH` variable

# Updating
- run `git pull` in the folder with code
- run build command for podman image (see Installation)
- make shell script executable

# Running
- Prepare the system using CHARMM-GUI. Generate GROMACS input files.
- Run the pipeline on gromacs folder from CHARMM-GUI: `bash run_chap.sh -i gromacs/ -o output`

# Troubleshooting
- Code might fail due to outdated GROMACS 2018. In particular, it fails to recognize `pcoupl` value `C-rescale`. To fix this, go inside GROMACS folder from CHARMM-GUI and run: `sed -i 's/C-rescale/Isotropic/g' *.mdp`. Instead of 'Isotropic' you can use: 'No' 'Berendsen' 'Parrinello-Rahman', 'MTTK’
- Program will overwrite any `mdp` file at destination. So edit input `mdp` files to adjust GROMACS behaviour.

# References
CHAP scripts for PyMOL are included with small adjustments