# Dockerfile optimized for matchairen-demo with HOMER
FROM continuumio/miniconda3:latest

# Install mamba in the base environment
RUN conda install -y mamba -n base -c conda-forge

# Install required Linux dependencies
RUN apt-get update && apt-get install -y \
    gawk \
    vim

# Copy HOMER to the container
COPY third_party/homer /opt/homer

# Add HOMER to PATH
ENV PATH="/opt/homer/bin:${PATH}"

# Copy the conda environment file to the container
COPY environment.yaml /tmp/environment.yaml

# Create conda environment using mamba (faster)
RUN mamba env create -f /tmp/environment.yaml

# Copy all pipeline code to the container
COPY . /work
WORKDIR /work

# Ensure scripts have execute permissions
RUN chmod +x run_snakemake.sh run_motif_pipeline.sh run_post_pipeline.sh prepare_data.bash

# Use bash and run the pipeline by default
SHELL ["/bin/bash", "-c"]
ENTRYPOINT ["bash", "-c", "source /opt/conda/etc/profile.d/conda.sh && conda activate GRN_software && ./run_snakemake.sh"]
