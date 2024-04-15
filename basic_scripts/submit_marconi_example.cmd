#!/bin/bash
#SBATCH --job-name=submit_example   # name of the job
#SBATCH --time=00:15:00             # walltime
#SBATCH --nodes=1                   # total number of nodes
#SBATCH --cpus-per-task=1           # Number of threads per task (=1, no OpenMP feature yet)
#SBATCH --ntasks-per-node=48        # MPI tasks (marconi, max 48)
#SBATCH --mem=64GB                  # memory requirement
#SBATCH --error=err.txt             # std-error file
#SBATCH --output=out.txt            # std-output file
#SBATCH --account=my_account        # account name
#SBATCH --partition=skl_fua_dbg     # partition name
#SBATCH --qos=normal                # quality of service

 # run the code
srun --cpu-bind=cores -l ./gyacomo.exe 2 24 1 0