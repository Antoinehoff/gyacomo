#!/bin/bash
#SBATCH --job-name=HeLaZ
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=64GB
#SBATCH --error=err.txt
#SBATCH --output=out.txt
#SBATCH --account=FUA36_TSVVT422
#SBATCH --partition=skl_fua_dbg
module load autoload hdf5 fftw
srun --cpu-bind=cores ./../../../bin/helaz3 2 12 2
