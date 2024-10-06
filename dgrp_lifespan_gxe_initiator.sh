#!/bin/bash
mkdir -p ./run/{log,logs_slurm} | sbatch ./dgrp_lifespan_gxe_snakemake_submitter.sh
