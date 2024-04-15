# This is an installation procedure for spcpc (Opensuze)
Clone the code

```bash
git clone https://gitlab.epfl.ch/ahoffman/gyacomo.git
```

Purge and load required modules
```bash
module purge
module load intel_comp/19.1 impi/5.0.3 mumps/4.10.0 hdf5/1.8.9
```

Enter `gyacomo` directory
```bash
cd gyacomo
```
make a lib directory somewhere (here we do it in `/gyacomo`)
```bash
mkdir lib
```
IMPORTANT: adapt the `gyacomo/local/dirs.inc` file by setting the path to your `gyacomo` folder in the `PREFIX` variable, and your path to the `lib` folder in the `LIBDIR` variable.

You can test now your installation and linkage by typing `make` in `/gyacomo`. The compilation should stop at the call of `futils` routines in `src/parallel_mod.F90`.

## Libraries
Install the following libraries in the `lib` folder.

### futils
Clone `futils` library (or ask for it if you do not have access)
```bash
git clone https://c4science.ch/diffusion/FUTILS/futils.git
```
go to `src`
```bash
cd /futils/src
```
adapt makefile
```bash
OPT = -O3 #-fast -Mvect=simd -Munroll -Mlarge_arrays
```
then make (still in `futils/src`)
```bash
make lib
```
put the lib files in directories
```bash
mkdir -p include/O
mkdir -p include/g
mkdir -p lib/O
mkdir -p lib/g
```
then install
```bash
make install
```
You can test now your installation and linkage by typing `make` in `/gyacomo`. The compilation should stop at the call of `fftw` routines in `src.fourier_mod.F90`.

### FFTW
Download `fftw` zip directory, unzip and go to the directory
```bash
wget http://www.fftw.org/fftw-3.3.10.tar.gz
tar -xvf fftw-3.3.10.tar.gz
cd fftw-3.3.10
```
Configure and make the double precision version
```bash
./configure --enable-mpi --prefix=$PWD
make
```
Install it (may require sudo command)
```bash
make install
```
Then same for the single precision
```bash
./configure --enable-float --enable-mpi --prefix=$PWD
make
```
Install it (may require sudo command)
```bash
make install
```
Gather all lib files in a new `lib` subdirectory
```bash
mkdir lib
cp lib64/* lib/.
cp *.la lib/.
```
The folder `/lib/fftw-3.3.10/lib` should contain
```bash
libfftw3.a    libfftw3f_mpi.a   libfftw3_mpi.a
libfftw3f.a   libfftw3f_mpi.la  libfftw3_mpi.la
libfftw3f.la  libfftw3.la
```
Go back to `/gyacomo` and verify that `FFTW3DIR   = $(LIBDIR)/fftw-3.3.10` in `/gyacomo/local/dirs.inc` when `ENVTYPE = Linux`.

You can test now your installation and linkage by typing `make` in `/gyacomo`. The compilation should stop at the call of `FM` routines in the `src/coeff_mod.F90` file.

### FM
Go back to your main `lib` directory (for us `/gyacomo/lib/`) and download the zipped folder there
```bash
wget https://dmsmith.lmu.build/FM1.4/FM_files.zip
```
Unzip and enter the folder
```bash
unzip FM_files.zip
mv FM_files FM
cd FM
```
First we have to copy the source file with a .F90 extensions (for our old intel compiler)
```bash
for file in *.f95 ; do cp "$file" "${file%.*}.F90" ; done
```
Now we compile manually (some of them are long to compile)
```bash
ifort fmsave.f95  -c -O3
ifort fm.f95  -c -O3
ifort fmzm90.f95  -c -O3
ifort TestFM.f95  -c -O3
ifort SampleFM.f95  -c -O3
ifort  fmsave.o  fm.o  fmzm90.o  TestFM.o  -o TestFM
ifort  fmsave.o  fm.o  fmzm90.o  SampleFM.o -o SampleFM
```
Then test the installation
```bash
./TestFM
./SampleFM
```
Put the library together and move it in a local `lib` directory
```bash
mkdir lib
ar r libfm.a fm.o fmsave.o fmzm90.o
mv libfm.a lib
```
Move the the `.mod` files in a `mod` local directory
```bash
mkdir mod
mv *.mod mod
```

You can test now your installation and linkage by typing `make` in `/gyacomo`. The compilation should go through and produce an executable in `/gyacomo/bin/.`named `gyacomo23_dp`.

## First runs
IF the executable `gyacomo23_dp` is present in `/gyacomo/bin/.`, you can now test a first run of the code.
In `/gyacomo`, run the simulations setup script
```bash
sh new_prob
```
This creates a new folder,`/gyacomo/simulations/` with an example of a simulation directory `/gyacomo/simulations/problem_01`. Inside, you can find a `tutorial.md` file explaining how to run the code and analyze its data.