FROM continuumio/miniconda3

# Avoid interactive apt prompts
ENV DEBIAN_FRONTEND=noninteractive

# --- System deps for ASCII video playback ---
RUN apt-get update && \
    apt-get install -y --no-install-recommends mpv libcaca0 ffmpeg locales && \
    rm -rf /var/lib/apt/lists/*

# (optional but nice) ensure a UTF-8 locale
RUN locale-gen en_US.UTF-8 || true
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

# --- Conda env ---
# build args for flexibility
ARG ENV_FILE=environment.yml
ARG ENV_NAME=bulk_rna_seq

# 1. Install Main Environment (Python + R + Tools)
#    Make sure 'seidr' is REMOVED from environment.yml before building!
COPY ${ENV_FILE} /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && conda clean -a

# 2. Install Seidr in an Isolated Environment
#    CHANGE: Use 'libgfortran=3' (no -ng) to get the legacy .so.3 library
RUN conda create -n seidr_env -c bioconda -c conda-forge -c defaults \
    seidr \
    "libgfortran=3" \
    && conda clean -a

# 3. Update PATH to include BOTH environments
#    Seidr binaries will be found here: /opt/conda/envs/seidr_env/bin
ENV PATH="/opt/conda/envs/${ENV_NAME}/bin:/opt/conda/envs/seidr_env/bin:${PATH}"

# --- MultiQC config (kept outside /app since it's a mount) ---
RUN mkdir -p /opt/multiqc
COPY multiqc_config.yaml /opt/multiqc/multiqc_config.yaml
ENV MULTIQC_CONFIG_PATH=/opt/multiqc/multiqc_config.yaml

# --- Source code ---
WORKDIR /usr/local/src/hulk
COPY . .

# install package in editable mode (Install into MAIN env)
RUN /opt/conda/envs/${ENV_NAME}/bin/pip install -e .

# --- Runtime workdir ---
WORKDIR /app

# put assets where the bind mount won't cover them
RUN mkdir -p /opt/hulk
COPY app/config/hulk_smash.mp4 /opt/hulk/hulk_smash.mp4

# entrypoint is the CLI installed by pip (points to app.cli:main)
ENTRYPOINT ["hulk"]
