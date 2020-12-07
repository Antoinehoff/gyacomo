#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10GB
#SBATCH --error=../results/Scaling/1024x512_L_100_Pe_1_Je_1_Pi_1_Ji_1_nB_0_nN_1_nu_1e-01_DG_mu_0e+00/err.txt
#SBATCH --output=../results/Scaling/1024x512_L_100_Pe_1_Je_1_Pi_1_Ji_1_nB_0_nN_1_nu_1e-01_DG_mu_0e+00/out.txt
#SBATCH --account=FUA34_GBSedge
#SBATCH --partition=skl_fua_dbg

#SBATCH --job-name=1024x512_L_100_Pe_1_Je_1_Pi_1_Ji_1_nB_0_nN_1_nu_1e-01_DG_mu_0e+00

module load intel
module load intelmpi
module load autoload hdf5/1.10.4--intelmpi--2018--binary
module load fftw
srun --cpu-bind=cores ./../bin/helaz