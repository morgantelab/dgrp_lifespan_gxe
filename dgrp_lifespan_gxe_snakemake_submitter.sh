#!/bin/bash

#SBATCH --job-name=dgrp_lifespan_gxe
#SBATCH --ntasks=1
#SBATCH --partition=compute
#SBATCH --time=336:00:00
#SBATCH --mem=2gb
#SBATCH --output=run/log/dgrp_lifespan_gxe.%j.out
#SBATCH --error=run/log/dgrp_lifespan_gxe.%j.err

source /data2/morgante_lab/fabiom/software/miniconda3/etc/profile.d/conda.sh
source /data2/morgante_lab/fabiom/software/miniconda3/etc/profile.d/mamba.sh
mamba activate dgrp_lifespan_gxe

#--dag | display | dot
#-p -n \
## test dag generation
#snakemake -p -n -s dgrp_lifespan_gxe_snakefile \
#          --configfile dgrp_lifespan_gxe.yaml \
#           --rerun-triggers mtime

snakemake \
  -s dgrp_lifespan_gxe_snakefile \
  --profile slurm \
  --latency-wait 120 \
  -k \
  --configfile dgrp_lifespan_gxe.yaml
