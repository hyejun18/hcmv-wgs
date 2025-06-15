#!/bin/zsh

#SBATCH --output=logs/alignment-%j.out
#SBATCH --mem-per-cpu=190G
source ~/.zshrc
mamba activate cgas
snakemake -s scripts/alignment.smk -j 128