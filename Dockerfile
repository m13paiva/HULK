FROM continuumio/miniconda3

# Suppress interactive prompts during apt-get package installations to ensure automated builds.
ENV DEBIAN_FRONTEND=noninteractive

# Install required system dependencies for ASCII video rendering and localization.
RUN apt-get update && \
    apt-get install -y --no-install-recommends mpv libcaca0 ffmpeg locales && \
    rm -rf /var/lib/apt/lists/*

# Generate and set UTF-8 locale to prevent character encoding errors during pipeline execution.
RUN locale-gen en_US.UTF-8 || true
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8
ENV MPLCONFIGDIR=/tmp

# ------------------------------------------------------------------------------
# Conda Environment Configuration
# ------------------------------------------------------------------------------
# Build arguments defining the environment specification file and target name.
ARG ENV_FILE=environment.yml
ARG ENV_NAME=bulk_rna_seq

# Initialize the primary environment (Python, R, and standard pipeline tools).
# Note: The 'seidr' package must be omitted from the environment.yml to prevent dependency conflicts.
COPY ${ENV_FILE} /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && conda clean -a

# Explicitly inject the BioConductor GO.db package into the primary environment.
RUN conda install -n ${ENV_NAME} -c bioconda -y bioconductor-go.db && conda clean -a

# Initialize an isolated Seidr environment.
# Seidr is installed in a separate environment using the legacy 'libgfortran=3' 
# to ensure compatibility with required .so.3 dynamic libraries.
RUN conda create -n seidr_env -c bioconda -c conda-forge -c defaults \
    seidr \
    "libgfortran=3" \
    && conda clean -a

# Expose binaries from both the primary environment and the isolated Seidr environment.
ENV PATH="/opt/conda/envs/${ENV_NAME}/bin:/opt/conda/envs/seidr_env/bin:${PATH}"

# ------------------------------------------------------------------------------
# Tool Configuration & Source Code Installation
# ------------------------------------------------------------------------------

# Deploy MultiQC configuration to an isolated path to prevent bind-mount overrides.
RUN mkdir -p /opt/multiqc
COPY multiqc_config.yaml /opt/multiqc/multiqc_config.yaml
ENV MULTIQC_CONFIG_PATH=/opt/multiqc/multiqc_config.yaml

# Transfer application source code to the container workspace.
WORKDIR /usr/local/src/hulk
COPY . .

# Install the primary pipeline package within the main Conda environment.
RUN /opt/conda/envs/${ENV_NAME}/bin/pip install -e .

# ------------------------------------------------------------------------------
# Runtime Preparation
# ------------------------------------------------------------------------------
WORKDIR /app

# Store static assets in a protected directory outside the standard runtime volume mounts.
RUN mkdir -p /opt/hulk
COPY app/config/hulk_smash.mp4 /opt/hulk/hulk_smash.mp4

# Define the default execution entrypoint mapped to the installed CLI application.
ENTRYPOINT ["hulk"]
