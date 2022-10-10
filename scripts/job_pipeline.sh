#! /bin/bash
#This script automatizes the launch of multiple job with chain dependency. 
nu_=0.01
dnu=0.01

Tm_=2000
dTm=2000

# First submit
lastji=$(sbatch submit_00.cmd)
lastjid=${lastji##* }
echo $lastji
#echo $lastjid

for i in {1..1}; do
    # Setup indices of job id (current and previous one)
    im1=$(awk "BEGIN {print $i-1}")
    idm1=$(printf '%02d' $im1)
    id=$(printf '%02d' $i)
    
    # Create new submit file from older one
    awk -v "ID=$id" '{ 
            if (NR == 8) print "#SBATCH --error=err_"ID".txt";
            else if (NR == 9) print "#SBATCH --output=out_"ID".txt";
            else if (NR == 12) print "srun --cpu-bind=cores ./gyacomo 2 24 1 "ID;
            else print $0}' submit_$idm1.cmd > submit_$id.cmd
 
    # Create new fort file from older one
    awk -v "NU=$nu_" -v "TM=$Tm_" -v "J2L=$im1" '{ 
            if (NR == 04) print "  tmax       = "TM;
       else if (NR == 40) print "  job2load   = "J2L;
       else if (NR == 54) print "  nu         = "NU;
       else print $0}' fort_$idm1.90 > fort_$id.90
    
    # Retrieve last jobid and launch next job with dep
    lastji=$(sbatch --dependency=afterok:$lastjid submit_0$i.cmd)
    lastjid=${lastji##* }
    echo $lastjid
    
    # Increment variables
    nu_=$(awk "BEGIN {print $nu_+$dnu}")
    Tm_=$(awk "BEGIN {print $Tm_+$dTm}")
done

