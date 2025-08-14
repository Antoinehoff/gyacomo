#! /bin/bash
#This script automatizes the launch of multiple job with chain dependency. 

# Function to display help
function show_help() {
    echo "Usage: $0 <simulation_id>"
    echo "This script submits the next job in the chain."
    echo "Ensure that mpi.config exists in the current directory."
    echo "The mpi.config file should contain 'Np Ny Nz'."
    echo "Options:"
    echo "  -h, --help         Show this help message"
    exit 1
}

# Check for help option
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
fi

# Check for correct usage
if [ "$#" -ne 1 ]; then
    echo "Error: Script requires exactly one argument."
    show_help
fi

if [ ! -f mpi.config ]; then
    echo "Error: mpi.config file not found."
    show_help
fi

# Read MPI configuration
read mpi1 mpi2 mpi3 < mpi.config
mpi_total=$((mpi1 * mpi2 * mpi3))

id=$1
idp1=$(awk "BEGIN {print $id + 1}")

# Format id and idp1 as two-digit integers
id_formatted=$(printf "%02d" $id)
idp1_formatted=$(printf "%02d" $idp1)

# Copy mpi.config to mpi_XX.config
cp mpi.config mpi_${idp1_formatted}.config

awk   -v "ID=$id_formatted" '{ 
        if (NR == 06) print "  job2load        = "ID;
else print $0}' params_${id_formatted}.in > params_${idp1_formatted}.in

awk   -v "IDP1=$idp1_formatted" -v "MPI_TOTAL=$mpi_total" -v "MPI1=$mpi1" -v "MPI2=$mpi2" -v "MPI3=$mpi3" '{ 
     if (NR == 07) print "#SBATCH --error=err_"IDP1".out";
     else if (NR == 08) print "#SBATCH --output=std_"IDP1".out";
     else if (NR == 12) print "srun -n "MPI_TOTAL" --cpu-bind=cores gyacomo23_dp "MPI1" "MPI2" "MPI3" " IDP1;
else print $0}' submit_${id_formatted}.sh > submit_${idp1_formatted}.sh

lastjid=$(cat jobid_${id_formatted}.txt)
echo sbatch --dependency=afterok:$lastjid submit_${idp1_formatted}.sh
submess=$(sbatch --dependency=afterok:$lastjid submit_${idp1_formatted}.sh)
jobid=${submess##* }
echo $jobid
echo $jobid > jobid_${idp1_formatted}.txt