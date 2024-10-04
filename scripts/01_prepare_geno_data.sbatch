#!/bin/bash

#SBATCH --job-name=prepare_geno_data
#SBATCH --cpus-per-task=1
#SBATCH --partition=fm-bigmem-1,fm-bigmem-2,fm-bigmem-3,fm-bigmem-4
#SBATCH --time=00:30:00
#SBATCH --mem=10G
#SBATCH --output=prepare_geno_data.%j.out
#SBATCH --error=prepare_geno_data.%j.err

###Activate environment
source /data2/morgante_lab/fabiom/software/miniconda3/etc/profile.d/conda.sh
source /data2/morgante_lab/fabiom/software/miniconda3/etc/profile.d/mamba.sh
mamba activate dgrp_lifespan_gxe

###Run command
plink --bfile /data2/morgante_lab/data/dgrp/genotypes/dgrp2 --maf 0.05 --geno 0.2 \
      --keep ../data/lines_kept.txt --make-bed --out ../data/dgrp_lifespan_gxe_filt_geno

###Deactivate environment
mamba deactivate
