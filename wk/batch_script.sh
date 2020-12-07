#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=32GB
#SBATCH --error=../results/Marconi/512x256_L_100_Pe_2_Je_1_Pi_2_Ji_1_nB_0.5_nN_1_nu_1e-01_FC_mu_5e-04/err.txt
#SBATCH --output=../results/Marconi/512x256_L_100_Pe_2_Je_1_Pi_2_Ji_1_nB_0.5_nN_1_nu_1e-01_FC_mu_5e-04/out.txt
#SBATCH --account=FUA34_GBSedge
#SBATCH --partition=skl_fua_prod

module load intel
module load intelmpi
module load autoload hdf5/1.10.4--intelmpi--2018--binary
module load fftw
srun --cpu-bind=cores ./../bin/helaz