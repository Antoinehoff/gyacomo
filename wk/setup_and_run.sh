#!/bin/bash
mkdir -p $CINECA_SCRATCH/HeLaZ/wk
mkdir -p $CINECA_SCRATCH/HeLaZ/bin

cd $CINECA_SCRATCH/HeLaZ/wk/
cp $HOME/HeLaZ/wk/fort.90 .
cp $HOME/HeLaZ/wk/batch_script.sh .
cp -r $HOME/HeLaZ/iCa ..

mkdir -p ../results/Marconi/512x256_L_100_Pe_2_Je_1_Pi_2_Ji_1_nB_0.5_nN_1_nu_1e-01_FC_mu_5e-04/
sbatch batch_script.sh
echo $CINECA_SCRATCH/HeLaZ/results/Marconi/512x256_L_100_Pe_2_Je_1_Pi_2_Ji_1_nB_0.5_nN_1_nu_1e-01_FC_mu_5e-04/out.txt