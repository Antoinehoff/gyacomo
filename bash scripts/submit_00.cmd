#!/bin/bash
#SBATCH --job-name=test_pipe
#SBATCH --time=00:02:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=64GB
#SBATCH --error=err_00.txt
#SBATCH --output=out_00.txt
#SBATCH --account=FUA36_TSVVT422
#SBATCH --partition=skl_fua_dbg
#SBATCH --qos=normal
srun --cpu-bind=cores ./gyacomo 2 24 1 0