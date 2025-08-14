#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

queue="regular"
account="m2116"
exe="gyacomo23_dp"
runtime="20:00:00"
mpi_decomposition="1 16 8"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -s)
            simid="$2"
            shift 2
            ;;
        -t)
            runtime="$2"
            shift 2
            ;;
        -m)
            mpi_decomposition="$2 $3 $4"
            shift 4
            ;;
        -q)
            queue="$2"
            shift 2
            ;;
        -a)
            account="$2"
            shift 2
            ;;
        -e)
            exe="$2"
            shift 2
            ;;
        -n)
            job_name="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 -s <simid> [-t <runtime>] [-m <np> <ny> <nz>] [-q <queue>] [-a <account>] [-e <executable>] [-n <job_name>]"
            echo "Options:"
            echo "  -s <simid>         Simulation ID (required)"
            echo "  -t <runtime>       Runtime in HH:MM:SS format (default: 20:00:00)"
            echo "  -m <np> <ny> <nz>  MPI decomposition (default: 1 32 8)"
            echo "  -q <queue>         Queue name (default: regular)"
            echo "  -a <account>       Account name (default: m2116)"
            echo "  -e <executable>    Executable name (default: gyacomo23_dp)"
            echo "  -n <job_name>      Custom job name (default: current directory name)"
            echo "  -h, --help         Display this help message and exit"
            exit 0
            ;;
        *)
            echo "Usage: $0 -s <simid> [-t <runtime>] [-m <np> <ny> <nz>] [-q <queue>] [-a <account>] [-e <executable>] [-n <job_name>]"
            exit 1
            ;;
    esac
done

check_job_status() {
    if [[ -z "$1" ]]; then
        echo "Usage: check_job_status <jobid>"
        return 1
    fi
    local jobid=$1
    local state
    state=$(squeue -j "$jobid" -h -o "%T" 2>/dev/null)
    if [[ -n "$state" ]]; then
        echo "$state"
        return 0
    fi
    state=$(sacct -j "$jobid" --format=State --noheader --parsable2 2>/dev/null | head -1)
    if [[ -n "$state" ]]; then
        echo "$state"
        return 0
    fi
    echo "NOT_FOUND"
    return 1
}

if [[ -z "$simid" ]]; then
    echo "Missing required argument: simid."
    echo "Usage: $0 -s <simid> [-t <runtime>] [-m <np> <ny> <nz>] [-q <queue>] [-a <account>] [-e <executable>] [-n <job_name>]"
    exit 1
fi

# Format simid as a two-digit number
simid=$(printf "%02d" $((10#$simid)))

# Use custom job name if provided, otherwise default to directory name
job_name=${job_name:-$(basename "$PWD")}

cores_per_node=$(lscpu | awk '/^Core\(s\) per socket:/ {cores=$NF} /^Socket\(s\):/ {sockets=$NF} END {print cores*sockets}')
err_file="err_${simid}.out"
out_file="std_${simid}.out"
job_file="jobid_${simid}.txt"
mpi_config="mpi_${simid}.config"

# Write MPI configuration to file
echo "$mpi_decomposition" > $mpi_config

read np ny nz <<< $(cat $mpi_config)
tasks=$((np * ny * nz))
nodes=$((tasks / cores_per_node))

slurm_script=$(cat <<EOF
#!/bin/bash
#SBATCH -N $nodes
#SBATCH -C cpu
#SBATCH -q $queue
#SBATCH -J $job_name
#SBATCH --account $account
#SBATCH --error=$err_file
#SBATCH --output=$out_file
#SBATCH -t $runtime
#SBATCH --cpus-per-task=1

srun -n $tasks --cpu-bind=cores $exe $np $ny $nz $simid
EOF
)

# Check for executable
if ! command -v "$exe" &> /dev/null && [[ ! -f "$exe" ]]; then
    echo "Error: Executable '$exe' not found in PATH or current directory."
    exit 1
fi

echo "$slurm_script"

read -p "Do you want to submit this job? [y]/n: " confirm
confirm=${confirm:-y}  # Default to 'y' if Enter is pressed
if [[ "$confirm" != "y" ]]; then
    echo "Job submission cancelled."
    exit 0
fi

# Check for params_XX.in file
simid_num=$((10#$simid))
prev_simid_num=$((simid_num - 1))
if [[ $simid_num -eq 0 ]]; then
    prev_simid_formatted="00"
else
    prev_simid_formatted=$(printf "%02d" $prev_simid_num)
fi
params_file="params_${prev_simid_formatted}.in"
new_params_file="params_${simid}.in"

if [[ ! -f "$params_file" ]]; then
    echo "Error: Required file $params_file not found."
    exit 1
fi

# Write the SLURM script to a file
slurm_script_file="submit_${simid}.sh"
echo "$slurm_script" > "$slurm_script_file"
chmod +x "$slurm_script_file"

if [[ $simid_num -gt 0 ]]; then
    # Copy and update job2load in the new param file
    awk -v "ID=$prev_simid_formatted" '{if ($1 == "job2load") print "  job2load   = "ID; else print $0}' "$params_file" > "$new_params_file"

    # Check for the previous job ID file
    prev_job_file="jobid_${prev_simid_formatted}.txt"
    if [[ ! -f "$prev_job_file" ]]; then
        echo "Error: Required file $prev_job_file not found."
        exit 1
    fi

    # Read the previous job ID
    prev_jobid=$(cat "$prev_job_file")
    prev_job_state=$(check_job_status "$prev_jobid")

    case "$prev_job_state" in
        "RUNNING"|"PENDING"|"COMPLETING")
            # Submit the job with a dependency
            job_output=$(sbatch --dependency=afterok:$prev_jobid "$slurm_script_file")
            dependency_info="Dependency: afterok:$prev_jobid"
            ;;
        *)
            # Submit the job without a dependency
            job_output=$(sbatch "$slurm_script_file")
            dependency_info=""
            ;;
    esac
else
    # Submit the job without a dependency
    job_output=$(sbatch "$slurm_script_file")
    dependency_info=""
fi

# Extract job ID and save it
jobid=$(echo "$job_output" | awk '{print $NF}')
echo "$jobid" > "$job_file"

# Print the SLURM submission message
echo "Job submitted successfully. SLURM output:"
echo "$job_output"
if [[ -n "$dependency_info" ]]; then
    echo "$dependency_info"
fi