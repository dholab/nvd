# Use Ubuntu 20.04 LTS as the base image
FROM ubuntu:20.04

# Set working directory
WORKDIR /scratch

# Set the maintainer label
LABEL maintainer="nrminor@wisc.edu"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York
ENV HOME=/opt
ENV ~=/opt

# Set default command to be the bash shell
ENTRYPOINT ["bash"]

# run a few apt installs
RUN apt-get update && \
    apt-get install -y curl wget git gcc g++ cmake util-linux && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir /dependencies && \
    dpkg -l > /dependencies/apt-get.lock

# Install all the dependencies locked with pyproject.toml:
# ----------------------------------------------------------------------------------- #
# 1) copy the required dependency and configuration file into the image
COPY pyproject.toml $HOME/pyproject.toml
COPY pixi.lock $HOME/pixi.lock

# 2) install pixi
RUN cd $HOME && PIXI_ARCH=x86_64 curl -fsSL https://pixi.sh/install.sh | bash

# 3) make sure pixi and pixi installs are on the $PATH
ENV PATH=$PATH:$HOME/.pixi/bin

# 4) install everything else with pixi
RUN cd $HOME && \
    script -q -c "pixi install --frozen" && \
    script -q -c "pixi clean cache --assume-yes" && \
    script -q -c "pixi add cxx-compiler cmake make"

# 5) modify the shell config so that each container launches within the pixi env
RUN echo "export PATH=$PATH:$HOME/.pixi/envs/default/bin" >> $HOME/.bashrc

# 6) modify some nextflow environment variables
RUN echo "export NXF_CACHE_DIR=/scratch" >> $HOME/.bashrc
RUN echo "export NXF_HOME=/scratch" >> $HOME/.bashrc

# ----------------------------------------------------------------------------------- #

# Copy necessary ncbi files
COPY conf/user-settings.mkfg /.ncbi/user-settings.mkfg

# Install NCBI tools and configure environment
RUN mkdir -p /.ncbi && \
    chmod 777 /.ncbi/user-settings.mkfg && \
    mkdir /build && cd /build && \
    git config --global http.sslVerify false && \
    git clone -b 3.2.0 https://github.com/ncbi/ngs-tools.git && \
    git clone -b 3.2.1 https://github.com/ncbi/ncbi-vdb.git && \
    git clone -b 3.2.1 https://github.com/ncbi/sra-tools.git && \
    sed -i 's/cmake_minimum_required[[:space:]]*([[:space:]]*VERSION[[:space:]]*2\.8\.12[[:space:]]*)/cmake_minimum_required(VERSION 3.5)/' /build/ngs-tools/CMakeLists.txt && \
    echo "Verifying CMake version change:" && \
    grep "cmake_minimum_required" /build/ngs-tools/CMakeLists.txt

# Build ncbi-vdb, sra-tools and ngs-tools
RUN cd /build/ncbi-vdb && \
    ./configure --relative-build-out-dir && \
    make -j$(nproc) && \
    cd /build/sra-tools && \
    ./configure --relative-build-out-dir && \
    make -j$(nproc) && \
    cd /build/ngs-tools && \
    ./configure --relative-build-out-dir && \
    make -j$(nproc) && \
    # Copy built binaries and clean up in the same layer
    find /build/OUTDIR -type f -executable -exec cp {} /usr/local/bin/ \; && \
    cd / && rm -rf /build
    
# remove now unnecessary compilers and clean the PyPI and conda caches
RUN cd $HOME && pixi remove rust cxx-compiler cmake make && pixi clean cache --yes

# Fix snakemake permission issue
RUN mkdir /.cache; chmod a+rwX /.cache


