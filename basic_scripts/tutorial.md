## Running Gyacomo Tutorial
This tutorial provides guidance on running Gyacomo in a new installation. It is designed to be used in conjunction with the `new_prob.sh` script, providing a refresher on the basics of Gyacomo commands.
It is meant to be followed once `sh new_prob.sh` has been run in `/gyacomo` which creates an example problem directory `/gyacomo/simulations/problem_01`. The following commands are expected to be run from the example problem directory.

### How to Run the Code Locally After Executing `new_prob.sh` in `/gyacomo` Directory?
- For a single-core run, execute:
    ```bash
    ./gyacomo.exe
    ```
- For a multi-core run, use:
    ```bash
    mpirun -np N ./gyacomo.exe np nky nz
    ```
    Here, `N` is the total number of processes, `np` is the number in the `p` direction (Hermite polynomials), `nky` is the number of poloidal wavenumbers, and `nz` is the number of parallel planes. The relation `N = np x nky x nz` must be satisfied.
In both commands above, the input parameters are expected to be in a `fort.90` file present in the directory where the code is run. You can append an additional integer `X` at the end of both single and multi-core runs, which points to an input file `fort_XX.90` where `XX` is the two-digit version of the number `X`. This allows for clearer restarts by indexing each consecutive run.

#### Examples
- Single-core run with parameters located in `fort_00.90`:
    ```bash
    ./gyacomo.exe 0
    ```
- Multi-core run with 2 cores in Hermite, 4 in ky, 1 in z, and reading the file `fort_00.90`:
    ```bash
    mpirun -np 8 ./gyacomo.exe 2 4 1 0
    ```
- Same as above but redirecting the standard terminal output to the file `out_00.txt` and running in the background:
    ```bash
    mpirun -np 8 ./gyacomo.exe 2 4 1 0 > out_00.txt &
    ```
    You can monitor the simulation progress with:
    ```bash
    tail -f out_00.txt
    ```

### How to Stop the Simulation?
Here are the stopping conditions of GYACOMO: 
- A nan (overflow) is detected in any field, this can be caused by a too high `dt` (CFL condition) or an unsufficient resolution (`Nkx`, `Nky`, `Nz`)
- The number of steps exceeds the number provided in `nrun`.
- The simulation run time exceeds the time defined in `maxruntime`.
- The simulation reaches the maximal physical time `tmax`.
- A file named `mystop`is present in the current directory.
Thus, you can smoothly stop your current simulation by creating an empty file named `mystop` in the directory where the code runs (here `gyacomo/simulations/problem_01`):
```bash
touch mystop
```
The code checks if such a file exists every 100 steps (this is set in the `gyacomo/src/tesend.F90` module). If it finds it, it will run the ending procedure and remove the file from the directory. 

### Checking the Output
Once the simulation is finished, the Python script `minimal_analysis.py` provides a basic example of result analysis. If you followed the tutorial, you should be able to obtain plots by typing:
```bash
python minimal_analysis.py
```

### Restarting a Simulation
In the parameter file, the parameter `job2load` defines from which previous outputs the code has to continue the simulation. If `job2load=-1`, a new start is made, and the output is located in the `outputs_00.h5` file in the current directory.

If the `outputs_XX` file exists, you can continue the run by setting `job2load=X` (where `XX` is the double-digit version of the integer `X`). It's recommended to link all runs with an input file `fort_XX.90`. In any case, you can find the input file used in the simulations in `outputs_XX.h5`, located in `outputs_XX.h5/files/STDIN.00`.

You can use `minimal_analysis.py X` to analyze the restart. The script will here look for the `outputs_XX.h5` file.

#### Example
This is a minimal example for a restart procedure
1. Copy the `fort_00.90` to a new `fort_01.90` input file.
2. Adapt the parameters of `fort_01.90` (Tmax, grads, anything).
3. **IMPORTANT:** Ensure that `fort_01.90` has `job2load = 0` (otherwise, it will restart a simulation from 0).
4. Now you can run:
    ```bash
    mpirun -np 8 ./gyacomo.exe 2 4 1 1 > out_01.txt &
    ```
    Note the `1` input in the fourth position; the code reads `fort_01.90`. We also direct the std output to a new file out_01.txt
5. You can monitor the run with:
    ```bash
    tail -f out_01.txt
    ```
5. Once the simulation is done (end or mystop), you can analyze the new data using:
    ```bash
    python minimal_analysis.py 1
    ```

### Running on Marconi

You can adapt and use the provided template `submit_marconi.cmd`.