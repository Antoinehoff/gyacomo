#! /bin/bash
#This script automatizes the launch of multiple job with possible chain dependency. 
SUBMIT=0   #to submit or not (use it to try first the file generation)
EXECPATH=gyacomo23_sp #path to the executable, will be renamed gyacomo in each subdir
Npp=1  #distribution of the processes (product must be 48 for one node usage)
Npy=12 #distribution of the processes (product must be 48 for one node usage)
Npz=4  #distribution of the processes (product must be 48 for one node usage)
Tm_=200    #maximal simulation time
runtime=10 #runtime asked to Marconi in hours
# program runtime in seconds -10% to have time to stop
RT_gyaco=$((3600 * runtime * 9 / 10))
# runtime to ask, in minutes
RT_slurm=$((60 * runtime))
# Number of time to continue the run (chain dep) N=0 means only one run
Ncont=0
# Hermite Laguerre resolution (the rest is taken from the fort_00.90 present in the current dir)
P=2
J=1
# Scan parameter names
P_1="delta"
P_2="kappa"
# Scan values
#A_1=(-0.3 -0.2 -0.1 0.0 0.1 0.2 0.3)
#A_2=(1.0 1.25 1.5 1.75 2.0)
A_1=(-0.3)
A_2=(1.0)
# find the lines where the parameters are found
L_1=$(grep -wn "$P_1" fort_00.90 | cut -d: -f1)
L_2=$(grep -wn "$P_2" fort_00.90 | cut -d: -f1)
# Check if the params exist
if [ ! -n "$L_1" ] || [ ! -n "$L_2" ] ; then
    echo "One of the parameter is not found in fort_00.90 -> quit"
    exit 1  # Exit the script with a non-zero status code
else
    echo "Scan over $P_1 (line $L_1) and $P_2 (line $L_2)"
fi
# Check the mpi distribution
Nptot=$((Npp * Npy * Npz))
# Check if the product is divisible by 48
if [ $((Nptot % 48)) -eq 0 ]; then
    Nnodes=$((Nptot / 48)) 
else
    echo "the distriubtion is not optimized for 48 processes"
    exit 1
fi
# loop over the indices
for i in "${!A_1[@]}"
do
	for j in "${!A_2[@]}" 
    do
        # Setup sim directory
        mkdir -p ${P_1}_${A_1[i]}_${P_2}_${A_2[j]}
        cd ${P_1}_${A_1[i]}_${P_2}_${A_2[j]}
        # Create new submit file from older one
        awk -v "P_=$P" -v "J_=$J" \
            -v "P1=${P_1}" -v "A1=${A_1[i]}" \
            -v "P2=${P_2}" -v "A2=${A_2[j]}" \
            -v "TR=$RT_slurm" -v "MPI=$MPIDISTR" \
            -v "np=$Npp" -v "ny=$Npy" -v "nz=$Npz" -v "Nnod=$Nnodes" '{
             if (NR == 2) print "#SBATCH --job-name=" P_ "_" J_ "_" P1 "_" A1 "_" P2 "_" A2; 
        else if (NR == 3) print "#SBATCH --time="TR;
        else if (NR == 4) print "#SBATCH --nodes="Nnod;
        else if (NR == 13)print "srun --cpu-bind=cores ./gyacomo " np " " ny " " nz " 0";
        else print $0}' ../submit_00.cmd > submit_00.cmd
        # Create new fort file from older one
        awk -v "P_=$P" -v "J_=$J" \
            -v "P1=${P_1}" -v "A1=${A_1[i]}" -v "L1=${L_1}" \
            -v "P2=${P_2}" -v "A2=${A_2[j]}" -v "L2=${L_2}" \
            -v "TM=$Tm_" -v "TR=$RT_gyaco" '{ 
             if (NR == 04) print "  tmax        = "TM;
        else if (NR == 05) print "  maxruntime  = "TR;
        else if (NR == 09) print "  pmax        = "P_;
        else if (NR == 10) print "  jmax        = "J_;
        else if (NR == L1) print "  "P1"        = "A1" ! Scanned";
        else if (NR == L2) print "  "P2"        = "A2" ! Scanned";
        else print $0}' ../fort_00.90 > fort_00.90
        cp ../$EXECPATH gyacomo
        if [ $SUBMIT -gt 0 ] ; then
            # Submit the job and display the msg
            submess=$(sbatch submit_00.cmd)
            jobid=${submess##* }
            echo $jobid
            # Save the jobid in jobid_00.txt
            echo $jobid > jobid_00.txt
            # Continue the job up to Ncontinue
            if [ $Ncont -gt 0 ]
            then
                a=0
                while [ $a -lt $Ncont ]
                do
                    sh ../continue.sh $a
                a=`expr $a + 1`
                done
            fi
        fi
        # Back in the previous directory
        cd ..
    done
done