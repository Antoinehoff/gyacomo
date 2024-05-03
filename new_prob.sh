#!/bin/bash
gyacomo_dir=$(pwd)

# Check if an argument is provided
if [ $# -ne 1 ]; then
    echo "- using 2D Z-pinch parameter"
    echo "  (you can ask for a cyclone base case by adding "CBC" after the call of the script)"
    inputfile=$gyacomo_dir/scripts/parameter_files/parameters_Zpinch.in
else
    # Check the value of the argument
    if [ "$1" = "help" ]; then
        echo "This script setup a new simulation directory in /simulation"
        echo "The basic setup is a small 2D z-pinch simulation."
        echo "You can ask a cyclone base case (more expensive) by typing sh new_prob.sh CBC"
        exit 1
    elif [ "$1" = "CBC" ]; then
        echo "- using CBC parameter"
        inputfile=$gyacomo_dir/scripts/parameter_files/parameters_CBC.in
    elif [ "$1" = "ZBC" ]; then
        echo "- using Zpinch parameter"
        inputfile=$gyacomo_dir/scripts/parameter_files/parameters_Zpinch.in
    else
        echo "ERROR STOP unknown entry"
        echo "(usage: sh new_prob.sh [CBC or ZBC (opt)])"
        exit 1
    fi
fi
echo "Setup of a new problem in the simulations directory..."

mkdir -p $gyacomo_dir/simulations
cd $gyacomo_dir/simulations

# Function to get the next available problem folder name
get_next_problem_folder() {
    local max_num=0
    for folder in problem_*; do
        num="${folder#problem_}"
        if [[ $num =~ ^[0-99]+$ ]] && [ "$num" -gt "$max_num" ]; then
            max_num="$num"
        fi
    done
    printf "problem_%02d" "$((max_num + 1))"
}

# Function to display usage information
usage() {
    echo "- usage: $0 [base_case_type] [folder_name]"
    echo "If folder_name is not provided, default is the next available 'problem_xx'"
    exit 1
}

# Check if an argument is provided
#if [ "$#" -gt 1 ]; then
#    usage
#fi

# Set the folder name based on the input or default
folder_name="${2:-$(get_next_problem_folder)}"

# Check if the executable exists
if [ ! -x $gyacomo_dir/bin/gyacomo23_dp ]; then
    echo "Error: Executable 'bin/gyacomo23_dp' not found or not executable. Please compile the code first."
    exit 1
fi

# Create the folder with the specified name if it doesn't exist
mkdir -p "$folder_name"

# Copy basic parameter file
cp $inputfile "$folder_name/params.in"

# Copy the marconi job submission example
cp $gyacomo_dir/scripts/parameter_files/submit_marconi_example.cmd "$folder_name/submit_marconi.cmd"

# Copy tutorial file
cp $gyacomo_dir/scripts/tutorials/tutorial.md "$folder_name/."

# Create a symbolic link to the executable
ln -s $gyacomo_dir/bin/gyacomo23_dp "$folder_name/gyacomo.exe"

# Create a symbolic link to the basic python analysis script
ln -s $gyacomo_dir/scripts/python_utilities/minimal_analysis.py "$folder_name/."
ln -s $gyacomo_dir/scripts/python_utilities/movie.py "$folder_name/."

echo "- folder 'simulations/$folder_name' created with symbolic link to the executable."
echo "- check the tutorial file in there!"
echo "...done"