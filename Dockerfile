# =============================================================================
# Dockerfile: GROMACS 2018 + GPU (CUDA) + CHAP
# Base: nvidia/cuda with Ubuntu 18.04
# =============================================================================


FROM nvidia/cuda:11.8.0-devel-ubuntu18.04


# NVIDIA OpenCL ICD is included with the CUDA image already.
# If using a plain Ubuntu base, install instead:
# nvidia-opencl-dev


LABEL maintainer="evgenii@gios.bio" \
      description="GROMACS 2018 with CUDA GPU support + CHAP (Channel annotation package)"

# ---------------------------------------------------------------------------
# 1. System dependencies
# ---------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        wget \
        curl \
        libfftw3-dev \
        libopenmpi-dev \
        openmpi-bin \
        libhdf5-dev \
        libboost-all-dev \
        libblas-dev \
        liblapack-dev \
        libgsl-dev \
        doxygen \
        python3 \
        python3-pip \
        python3-dev \
        pkg-config \
        zlib1g-dev \
        libssl-dev \
        ca-certificates \
        ocl-icd-opencl-dev \
        ocl-icd-libopencl1 \
        opencl-headers \
        clinfo \
        libboost-all-dev \
        libblas-dev liblapacke-dev libatlas-base-dev libopenblas-dev \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# 2. Build & install GROMACS 2018.8 (latest 2018 patch) with CUDA support
# --------------------------------------------------------------------------
ENV GROMACS_VERSION=2018
ENV GROMACS_INSTALL_DIR=/opt/gromacs

WORKDIR /tmp

RUN wget -q http://ftp.gromacs.org/pub/gromacs/gromacs-${GROMACS_VERSION}.tar.gz \
    && tar xzf gromacs-${GROMACS_VERSION}.tar.gz \
    && mkdir -p gromacs-${GROMACS_VERSION}/build \
    && cd gromacs-${GROMACS_VERSION}/build \
    && cmake .. \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=${GROMACS_INSTALL_DIR} \
        -DGMX_GPU=OpenCL \
        -DGMX_USE_OPENCL=ON \
        -DGMX_FFT_LIBRARY=fftw3 \
        -DGMX_BUILD_OWN_FFTW=OFF \
        -DREGRESSIONTEST_DOWNLOAD=OFF \
        -DGMX_DOUBLE=OFF \
    && make -j$(nproc) \
    && make install \
    && cd /tmp \
    && rm -rf gromacs-${GROMACS_VERSION}*

# Source GROMACS environment for all subsequent commands
ENV PATH="${GROMACS_INSTALL_DIR}/bin:${PATH}" \
    LD_LIBRARY_PATH="${GROMACS_INSTALL_DIR}/lib:${LD_LIBRARY_PATH}"

# ---------------------------------------------------------------------------
# 3. Build & install CHAP (Channel Annotation Package)
#    
# ---------------------------------------------------------------------------
ENV CHAP_VERSION=0_9_1
ENV CHAP_INSTALL_DIR=/opt/chap

RUN wget https://github.com/channotation/chap/archive/version_${CHAP_VERSION}.tar.gz \
    && mkdir -p /tmp/chap/build \
    && tar xfvz version_${CHAP_VERSION}.tar.gz -C chap --strip-components 1 \
    && cd /tmp/chap/build \
    && cmake .. \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_PREFIX_PATH=${GROMACS_INSTALL_DIR} \
        -DCMAKE_INSTALL_PREFIX=${CHAP_INSTALL_DIR} \
        -DGROMACS_DIR=${GROMACS_INSTALL_DIR} \
    && make -j$(nproc) \
    && make install \
    && rm -rf /tmp/chap

ENV PATH="${CHAP_INSTALL_DIR}/bin:${PATH}" \
    LD_LIBRARY_PATH="${CHAP_INSTALL_DIR}/lib:${LD_LIBRARY_PATH}"

# ---------------------------------------------------------------------------
# 4. Shell initialisation – source GROMACS env on login
# ---------------------------------------------------------------------------
RUN echo "source ${GROMACS_INSTALL_DIR}/bin/GMXRC" >> /etc/bash.bashrc

# ---------------------------------------------------------------------------
# 5. Working directory
# ---------------------------------------------------------------------------
WORKDIR /workspace

CMD ["/bin/bash"]
