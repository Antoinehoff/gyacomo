#! /bin/bash
#This script automatizes the launch of multiple job with chain dependency. 
id=$1
idp1=$(awk "BEGIN {print $id + 1}")
    awk   -v "ID=$id" '{ 
            if (NR == 06) print "  job2load        = "ID;
    else print $0}' fort_0$id.90 > fort_0$idp1.90

    awk   -v "IDP1=$idp1" '{ 
         if (NR == 08) print "#SBATCH --error=err_0"IDP1".txt";
         else if (NR == 09) print "#SBATCH --output=out_0"IDP1".txt";
         else if (NR == 12) print "srun --cpu-bind=cores ./../gyacomo23_sp 8 6 4 " IDP1
    else print $0}' submit_0$id.cmd > submit_0$idp1.cmd

    lastjid=$(cat jobid_0$id.txt)
    #echo sbatch --dependency=afterok:$lastjid submit_0$idp1.cmd
    #submess=$(sbatch --dependency=afterok:$lastjid submit_0$idp1.cmd)
    jobid=${submess##* }
    echo $jobid
    echo $jobid > jobid_0$idp1.txt