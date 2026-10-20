# Helper image for processes that run workflow/scripts/*.py (pandas/openpyxl) and FETCH_SRA (sra-tools).
# Build:  docker build -t ghcr.io/<owner>/blactx-helper:0.1.0 containers/helper
FROM mambaorg/micromamba:2.0.5
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.12 pandas openpyxl biopython sra-tools procps-ng \
    && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
