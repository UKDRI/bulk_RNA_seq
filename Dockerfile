FROM continuumio/miniconda3:latest
WORKDIR /app
RUN conda install -y python=3.11 && conda clean -a -y
RUN conda update -n base -c defaults conda -y && \
    conda install -y -c conda-forge -c bioconda snakemake=8.9.0 graphviz && \
    conda clean -a -y
COPY . /app
ENTRYPOINT ["snakemake"]
CMD ["--help"]
