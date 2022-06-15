# Set the base image to Ubuntu
FROM continuumio/miniconda3:latest

# Install tools using conda
RUN conda config --add channels bioconda && conda config --add channels conda-forge && \
  conda update -n base -c defaults --yes conda && apt-get update && apt-get -y install gcc make
  
RUN conda install mamba 

RUN mamba install --yes hisat2 gatk4 "samtools>=1.10" trim-galore bedtools bamtools bwa bbmap "openssl>=1.0" samblaster bcftools rtg-tools

RUN mamba update -y -c r ncurses 
