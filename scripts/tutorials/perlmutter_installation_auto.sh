#! /bin/bash
# This tutorial helps you install Gyacomo on Perlmutter.

# Record the start time
start=$(date +%s)

# Display title
echo '=============================================='
echo '   ______                                     '
echo '  / ____/_  ______ __________  ____ ___  ____ '
echo ' / / __/ / / / __ `/ ___/ __ \/ __ `__ \/ __ \'
echo '/ /_/ / /_/ / /_/ / /__/ /_/ / / / / / / /_/ /'
echo '\____/\__, /\__,_/\___/\____/_/ /_/ /_/\____/ '
echo '     /____/                                   '
echo '---This is the GYACOMO code auto installer----'
echo 'for Perlmutter (local HDF5, futils, FFTW)'
echo 'takes ~5min                                   '
echo '=============================================='

echo "1/6 install HDF5"

echo "2/6 install futils"

echo "3/6 install FM"

echo "2/6 install fftw"



# Calculate execution time
execution_time=$((end - start))

# Check if execution time is more than a minute
if [ $execution_time -gt 60 ]; then
    execution_time_minutes=$((execution_time / 60))
    execution_time_seconds=$((execution_time % 60))
    echo "Script executed in $execution_time_minutes minute(s) and $execution_time_seconds second(s)."
else
    echo "Script executed in $execution_time second(s)."
fi
echo '=============================================='