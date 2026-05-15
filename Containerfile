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

# run a few apt installs
RUN apt-get update && \
    apt-get install -y curl util-linux && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir /dependencies && \
    dpkg -l > /dependencies/apt-get.lock

# Install all the dependencies locked with pyproject.toml:
# ----------------------------------------------------------------------------------- #
# 1) copy the required dependency and configuration files into the image
COPY pyproject.toml $HOME/pyproject.toml
COPY pixi.lock $HOME/pixi.lock
COPY uv.lock $HOME/uv.lock
COPY lib/ $HOME/lib/
COPY main.nf $HOME/main.nf
COPY nextflow.config $HOME/nextflow.config

# 2) install pixi
RUN cd $HOME && PIXI_ARCH=x86_64 curl -fsSL https://pixi.sh/install.sh | bash

# 3) make sure pixi and pixi installs are on the $PATH
ENV PATH=$PATH:$HOME/.pixi/bin

# 4) install everything else with pixi
RUN cd $HOME && \
    script -q -c "pixi install --frozen" && \
    script -q -c "pixi clean cache --assume-yes"

# 5) Add pixi environment to PATH (works in Docker, Podman, AND Apptainer)
ENV PATH=$PATH:/opt/.pixi/envs/default/bin

# 6) Set Nextflow environment variables (works in all container runtimes)
ENV NXF_CACHE_DIR=/scratch
ENV NXF_HOME=/scratch

# Fix snakemake permission issue
RUN mkdir /.cache; chmod a+rwX /.cache
