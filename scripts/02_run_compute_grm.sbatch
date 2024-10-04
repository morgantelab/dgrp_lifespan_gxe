#!/bin/bash

#SBATCH --job-name=run_compute_grm
#SBATCH --cpus-per-task=4
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=00:30:00
#SBATCH --mem=15G
#SBATCH --output=run_compute_grm.%j.out
#SBATCH --error=run_compute_grm.%j.err

###Load module
source /data2/morgante_lab/fabiom/software/miniconda3/etc/profile.d/conda.sh
source /data2/morgante_lab/fabiom/software/miniconda3/etc/profile.d/mamba.sh
mamba activate dgrp_lifespan_gxe

###Set variables
export MKL_NUM_THREADS=4
export OMP_NUM_THREADS=4

###Run script
Rscript ../code/compute_grm.R --geno ../data/dgrp_lifespan_gxe_filt_geno \
                              --output ../data/dgrp_lifespan_gxe_grm.rds

###Unload module
mamba deactivate
