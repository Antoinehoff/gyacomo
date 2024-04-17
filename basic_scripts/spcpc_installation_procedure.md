This tutorial helps you install Gyacomo on an SPC medusa PC, based on OpenSUSE and intel compiler. You will be guided from the cloning of the repository to the first run.

## 1. Prelude
1.1 - First, we clone the repository.
```bash
git clone https://gitlab.epfl.ch/ahoffman/gyacomo.git
```

1.2 - Then we purge and load the following modules (the default one may be of older intel compiler version, which will not work with FFTW3).
```bash
module purge
module load intel_comp/19.1 impi/5.0.3 mumps/4.10.0 hdf5/1.8.9
```

1.3 - Now we navigate to the `gyacomo` directory.
```bash
cd gyacomo
```

1.4 - Finally, we make a `lib` directory (here we make it in `/gyacomo` but it can be anywhere as long as the path to it is set correctly in `gyacomo/local/dirs.inc`).
```bash
mkdir lib
```
**Important:** Verify that the `gyacomo/local/dirs.inc` file points to your `gyacomo` folder in the `PREFIX` variable and to your `lib` folder in the `LIBDIR` variable.

1.5 - We can do a first test the installation and linkage by typing `make` in `/gyacomo`. The compilation should halt with an error at the call of `futils` routines in `src/parallel_mod.F90`. This is normal, let us continue.

## Install Libraries
Now you will be guided to install the three libraries required by Gyacomo:
- Futils, a wrapper library to handle easily h5 outputs (developped by T.M.Tran at SPC)
- FM, an arbitrary precision library used to compute the Laguerre-Laguerre product coefficients `$d_{njs}$` without accuracy loss (see [coeff module](Sources/coeff)).
- FFTW3, the Fastest Fourier Transform in the West version 3, which enables us to compute the nonlinear term using a pseudo-spectral approach (see [fourier module](Sources/fourier)).
Each of these libraries must be installed in the `/lib` folder (here `/gyacomo/lib`).

### 2. Futils installation
2.1 - We navigate to our `lib` folder and clone first the `futils` library.
```bash
cd lib
git clone https://c4science.ch/diffusion/FUTILS/futils.git
```
**Note:** if you do not have the access, contact me and I will provide a .zip file.

2.2 - Then, we navigate to the `src` directory.
```bash
cd futils/src
```

2.3 - You have to adapt the makefile because some compilation options are not available with the intel compiler (we put a simple -O3 option).
```bash
sed -i '37s/.*/OPT = -O3/' Makefile
```
**Note:** the Makefile is also assuming that the Hdf5 library path is located in a `$(HDF5)` variable and that the `mpif90` compiler is defined. If you are on Marconi HPCC, you may have to change line 34 to `F90 = mpiifort` and line 35 to `HDF5 = $(HDF5_HOME)` (having prealably loaded the Hdf5 module)

2.4 - We can now compile the library.
```bash
make lib
```

2.5 - We create necessary directories to store our library files.
```bash
mkdir -p include/O
mkdir -p include/g
mkdir -p lib/O
mkdir -p lib/g
```

2.6 - Time now to install the library (this will put the library files in the previously made directories).
```bash
make install
```

2.7 - We can test now the installation and linkage of futils by typing `make` in `/gyacomo`. The compilation should pass the previously observed error and halt now with an error at the call of `fftw` routines in `src.fourier_mod.F90`.

### 3. FFTW3 installation
3.1 - We download the `fftw` zip directory, unzip it, and navigate to the directory.
```bash
wget http://www.fftw.org/fftw-3.3.10.tar.gz
tar -xvf fftw-3.3.10.tar.gz
cd fftw-3.3.10
```

3.2 - We configure the installation and make the double-precision version.
```bash
./configure --enable-mpi --prefix=$PWD
make
```

3.3 -  We can now install it.
```bash
make install
```

3.4 - For single-precision runs, we also configure and install the single-precision version.
```bash
./configure --enable-float --enable-mpi --prefix=$PWD
make
make install
```

3.5 - We gather the library files in a new `lib` subdirectory.
```bash
mkdir lib
cp lib64/* lib/.
cp *.la lib/.
```

3.7 - We go back to `/gyacomo` and verify that `FFTW3DIR   = $(LIBDIR)/fftw-3.3.10` in `/gyacomo/local/dirs.inc` when `ENVTYPE = Linux`.

3.8 - We test the installation and linkage by typing `make` in `/gyacomo`. The compilation should halt with an error at the call of `FM` routines in the `src/coeff_mod.F90` file.

## 4. FM installation
4.1 - Navigate back to the main `lib` directory (e.g., `/gyacomo/lib/`) and download the zipped folder.
```bash
wget https://dmsmith.lmu.build/FM1.4/FM_files.zip
```

4.2 - Unzip and enter the folder.
```bash
unzip FM_files.zip
mv FM_files FM
cd FM
```

4.3 - Copy the source file with a `.f95` to `.F90`extension. This is required to use the intel compiler (gfortran compilation here will not work if Gyacomo is compiled with mpiifort).
```bash
for file in *.f95 ; do cp "$file" "${file%.*}.F90" ; done
```

4.4 - We compile manually (some compilation may take ~5min).
```bash
ifort fmsave.f95 -c -O3
ifort fm.f95 -c -O3
ifort fmzm90.f95 -c -O3
ifort TestFM.f95 -c -O3
ifort SampleFM.f95 -c -O3
ifort fmsave.o fm.o fmzm90.o TestFM.o -o TestFM
ifort fmsave.o fm.o fmzm90.o SampleFM.o -o SampleFM
```

4.5 We can test the installation.
```bash
./TestFM
./SampleFM
```

4.6 The library file is put together and moved to a local `lib` directory.
```bash
mkdir lib
ar r libfm.a fm.o fmsave.o fmzm90.o
mv libfm.a lib
```

4.7 Move the `.mod` files to a `mod` local directory.
```bash
mkdir mod
mv *.mod mod
```

4.8 - We test the installation and linkage by typing `make` in `/gyacomo`. The compilation should proceed and produce an executable in `/gyacomo/bin/.` named `gyacomo23_dp`.

## 5. First Run
5.1 - The executable `gyacomo23_dp` should be present in `/gyacomo/bin/.`. We can now test the first run of the code. In `/gyacomo`, run the simulations setup script:
```bash
sh new_prob.sh ZBC
```
`ZBC`stands for Z-pinch base case, which is a 2D nonlinear simulation (recommended for a light run). You can also call the script with `CBC`, which is a minimal example of a cyclone base case (underresolved to be light).

This creates a new folder, `/gyacomo/simulations/`, with an example of a simulation directory `/gyacomo/simulations/problem_01`. Here we present a basic method to run it. You can find a `tutorial.md` file explaining how to run the code and analyze its data in the simulation folder.

5.1 - Go to the simulation folde `/gyacomo/simulations/problem_01/' and run the code
```bash
cd /gyacomo/simulations/problem_01/
mpirun -np 4 ./gyacomo 1 4 1 0 > out_00.txt &
```

5.2 - The code is running now in background. You can check the advancement typing
```bash
tail -f out_00.txt
```
(`ctrl+c` stop following the std output`)

5.3 - The run will stop once Tmax is reached in the fort_00.90 parameter file. However, you can stop it earlier by typing
```bash
touch mystop
````
which creates an empty file in the simulation folder that acts as a stopping signal for the code (you can also use `pkill gyacomo` but this will make the output data unusable).

5.4 - Once the run is finished, or stoped with `mystop`, you can analyze the data with a minimal python analysis script.
```bash
python minimal_analysis.py
```
This will provide the time traces of heat and particle fluxes as well as the last frame of the ion gyrodensity and electrostatic potential in the outboard midplane.
