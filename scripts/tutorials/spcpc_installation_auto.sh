#! /bin/bash
# This tutorial helps you install Gyacomo on an SPC medusa PC, based on OpenSUSE and intel compiler. You will be guided from the cloning of the repository to the first run.

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
echo 'This is the GYACOMO code auto installer'
echo '=============================================='

# 1. Prelude
# 1.1 - First, we clone the repository.
echo "1/6 Clone the directory "
git clone https://gitlab.epfl.ch/ahoffman/gyacomo.git &> gyacomo_clone.out
cd gyacomo

# 1.2 - Then we purge and load the following modules (the default one may be of older intel compiler version, which will not work with FFTW3).
echo "      Purge and load the modules"
module purge
module load intel_comp/19.1 impi/5.0.3 mumps/4.10.0 hdf5/1.8.9

# 1.3 - Now we navigate to the `gyacomo` directory.
# cd gyacomo

# 1.4 - Finally, we make a `lib` directory (here we make it in `/gyacomo` but it can be anywhere as long as the path to it is set correctly in `gyacomo/local/dirs.inc`).
echo "      make lib directory"
mkdir lib
# **Important:** Verify that the `gyacomo/local/dirs.inc` file points to your `gyacomo` folder in the `PREFIX` variable and to your `lib` folder in the `LIBDIR` variable.

# 1.5 - We can do a first test the installation and linkage by typing `make` in `/gyacomo`. The compilation should halt with an error at the call of `futils` routines in `src/parallel_mod.F90`. This is normal, let us continue.

## Install Libraries
# Now you will be guided to install the three libraries required by Gyacomo:
# - Futils, a wrapper library to handle easily h5 outputs (developped by T.M.Tran at SPC)
# - FM, an arbitrary precision library used to compute the Laguerre-Laguerre product coefficients `$d_{njs}$` without accuracy loss (see [coeff module](Sources/coeff)).
# - FFTW3, the Fastest Fourier Transform in the West version 3, which enables us to compute the nonlinear term using a pseudo-spectral approach (see [fourier module](Sources/fourier)).
# Each of these libraries must be installed in the `/lib` folder (here `/gyacomo/lib`).

echo "      repository cloned"
echo "______________________________________________"

#############################################
### 2. Futils installation
# 2.1 - We navigate to our `lib` folder and clone first the `futils` library.
echo "2/6 Installing Futils"
cd lib
echo "      clone repository"
git clone https://c4science.ch/diffusion/FUTILS/futils.git &> futils_clone.out
# **Note:** if you do not have the access, contact me and I will provide a .zip file.

# 2.2 - Then, we navigate to the `src` directory.
cd futils/src

# 2.3 - You have to adapt the makefile because some compilation options are not available with the intel compiler (we put a simple -O3 option).
sed -i '37s/.*/OPT = -O3/' Makefile
# **Note:** the Makefile is also assuming that the Hdf5 library path is stored in a `$(HDF5)` variable and that the `mpif90` compiler is defined. If you are on Marconi HPCC, you may have to change line 34 to `F90 = mpiifort` and line 35 to `HDF5 = $(HDF5_HOME)` (having prealably loaded the Hdf5 module)

# 2.4 - We can now compile the library.
echo "      make library.."
make lib &> futils_make_lib.out

# 2.5 - We create necessary directories to store our library files.
mkdir -p include/O
mkdir -p include/g
mkdir -p lib/O
mkdir -p lib/g

# 2.6 - Time now to install the library (this will put the library files in the previously made directories).
echo "      compiling.."
make install &> futils_make_install.out

# 2.7 - We can test now the installation and linkage of futils by typing `make` in `/gyacomo`. The compilation should pass the previously observed error and halt now with an error at the call of `fftw` routines in `src.fourier_mod.F90`.
cd ../../
echo "      Futils installed"
echo "______________________________________________"

#############################################
### 3. FFTW3 installation
echo "3/6 Installing FFTW3"
echo "      download source"
# 3.1 - We download the `fftw` zip directory, unzip it, and navigate to the directory.
wget http://www.fftw.org/fftw-3.3.10.tar.gz &> wget_fftw3.out
tar -xvf fftw-3.3.10.tar.gz &> tar_fftw3.out
cd fftw-3.3.10

# 3.2 - We configure the installation and make the double-precision version.
echo "      configure double precision mpi.."
./configure --enable-mpi --prefix=$PWD &> configure_mpi.out
echo "      compile double precision mpi.."
make &> make.out

# 3.3 -  We can now install it.
echo "      installing.."
make install &> make_install.out

# 3.4 - For single-precision runs, we also configure and install the single-precision version.
echo "      configure single precision mpi.."
./configure --enable-float --enable-mpi --prefix=$PWD &> configure_float.out
echo "      compile single precision mpi.."
make &> make_float.out
echo "      installing.."
make install &> make_install_float.out

# 3.5 - We gather the library files in a new `lib` subdirectory.
mkdir lib
cp -r lib64/* lib/.
cp -r *.la lib/.

# 3.7 - We go back to `/gyacomo` and verify that `FFTW3DIR   = $(LIBDIR)/fftw-3.3.10` in `/gyacomo/local/dirs.inc` when `ENVTYPE = Linux`.
cd ../

# 3.8 - We test the installation and linkage by typing `make` in `/gyacomo`. The compilation should halt with an error at the call of `FM` routines in the `src/coeff_mod.F90` file.
echo "      FFTW3 installed"
echo "______________________________________________"


#############################################
## 4. FM installation
echo "4/6 Installing FM"

# 4.1 - Navigate back to the main `lib` directory (e.g., `/gyacomo/lib/`) and download the zipped folder.
echo "      download source"
wget https://dmsmith.lmu.build/FM1.4/FM_files.zip &> wget_FM.out

# 4.2 - Unzip and enter the folder.
unzip FM_files.zip &> unzip_FM.out
mv FM_files FM
cd FM

# 4.3 - Copy the source file with a `.f95` to `.F90`extension. This is required to use the intel compiler (gfortran compilation here will not work if Gyacomo is compiled with mpiifort).
for file in *.f95 ; do cp "$file" "${file%.*}.F90" ; done

# 4.4 - We compile manually (some compilation may take ~5min).
echo "      compiling.. (~5min)"
echo "      fmsave.."
ifort fmsave.F90 -c -O3
echo "      fm.."
ifort fm.F90 -c -O3
echo "      fmzm.."
ifort fmzm90.F90 -c -O3
# echo "      TestFM.."
# ifort TestFM.F90 -c -O3
# echo "      SampleFM.."
# ifort SampleFM.F90 -c -O3
# echo "      building all.."
# ifort fmsave.o fm.o fmzm90.o TestFM.o -o TestFM
# ifort fmsave.o fm.o fmzm90.o SampleFM.o -o SampleFM
echo "      ..done"
# 4.5 We can test the installation.
# ./TestFM &> TestFM.out
# ./SampleFM &> SampleFM.out

# 4.6 The library file is put together and moved to a local `lib` directory.
mkdir lib
ar r libfm.a fm.o fmsave.o fmzm90.o
mv libfm.a lib

# 4.7 Move the `.mod` files to a `mod` local directory.
mkdir mod
mv *.mod mod

cd ../
echo "      FM installed"
echo "______________________________________________"

#############################################
# 5 compilation of GYACOMO
cd ../
echo "5/6 Compilation of GYACOMO"
make 2>&1 | tee make.out
echo "______________________________________________"

#############################################
# 6. setup of a testcase
echo "6/6 Test case setup"
make new_prob
echo "-- GYACOMO is ready --"
echo "you can test running it using:"
echo "cd gyacomo/simulations/problem_01"
echo "mpirun -np 4 ./gyacomo.exe 1 4 1"
echo "______________________________________________"

#############################################
# 7. Epilogue
echo "End of the script"
# Record end time
end=$(date +%s)
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