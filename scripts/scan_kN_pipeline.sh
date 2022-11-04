#! /bin/bash
# lastjid=$(sbatch submit_00.cmd)

nu_=0.01
dnu=0.0
kn_=1.7
dkn=0.2
kt_=0.425
dkt=0.05
Lx_=150
dLx=030

Tm_=4000
dTm=1000

# First submit
istart=0
lastji=$(sbatch submit_0$istart.cmd)
lastjid=${lastji##* }
echo $lastji

for i in {1..5}; do
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
    awk -v "NU=$nu_" -v "TM=$Tm_" -v "J2L=$im1" -v "KN=$kn_" -v "KT=$kt_" -v "LX=$Lx_" '{ 
            if (NR == 04) print "  tmax       = "TM;
       else if (NR == 40) print "  job2load   = "J2L;
       else if (NR == 13) print "  Lx         = "LX;
       else if (NR == 54) print "  nu         = "NU;
       else if (NR == 61) print "  K_Ne       = "KN;
       else if (NR == 62) print "  K_Ni       = "KN;
       else if (NR == 63) print "  K_Te       = "KT;
       else if (NR == 64) print "  K_Ti       = "KT;
       else print $0}' fort_$idm1.90 > fort_$id.90
    
    # Retrieve last jobid and launch next job with dep
    lastji=$(sbatch --dependency=afterok:$lastjid submit_$id.cmd)
    lastjid=${lastji##* }
    echo $lastjid
    
    # Increment variables
    Lx_=$(awk "BEGIN {print $Lx_+$dLx}")
    nu_=$(awk "BEGIN {print $nu_+$dnu}")
    kn_=$(awk "BEGIN {print $kn_+$dkn}")    
    kt_=$(awk "BEGIN {print $kt_+$dkt}")
    Tm_=$(awk "BEGIN {print $Tm_+$dTm}")
done

