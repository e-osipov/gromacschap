# gromacschap
Gromacs and CHAP pipeline  automated pipeline

# Installation
- Install Podman >4.2 wtih GPU support (ask Claude for specific commands for your distro)
- fetch code using `git clone`
- Adjust `Dockerfile` for your architecture:
    - first, find out what is your GPU architecture. Code for Nvidia GPU only: `nvidia-smi --query-gpu=gpu_name,compute_cap --format=csv`. For example, RTX4000 has compute 8.9
    - Multipy number from previous step by 10 and change line in `Dockerfile`: `-DGMX_CUDA_TARGET_COMPUTE="89" \` to show correct value
- Build image: `podman build  --format docker -t chap:1.0 -t chap:latest .`
- Once done, run the pipeline:
