#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -J gyacomo_test
#SBATCH --account m2116
#SBATCH --output std.out
#SBATCH --error err.out
#SBATCH -t 00:30:00

#run the application:
srun -n 32 --cpu_bind=cores gyacomo.exe 1 32 1
exit 0