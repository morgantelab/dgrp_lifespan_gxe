#!/bin/bash

#SBATCH --job-name=run_plot_varcomps_gencor
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --output=run_plot_varcomps_gencor.%j.out
#SBATCH --error=run_plot_varcomps_gencor.%j.err

###Load module

source /data2/morgante_lab/fabiom/software/miniconda3/etc/profile.d/conda.sh
source /data2/morgante_lab/fabiom/software/miniconda3/etc/profile.d/mamba.sh
mamba activate dgrp_lifespan_gxe

###Set variables
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Run script
Rscript ../code/paper_figures/plot_varcomps_gencor.R

###Unload module
mamba deactivate
