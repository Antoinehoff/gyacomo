# Installation Procedure for SPCPC (OpenSUSE)

## 1. Prelude
1.1 Clone the repository:
```bash
git clone https://gitlab.epfl.ch/ahoffman/gyacomo.git
```

1.2 Purge and Load Required Modules
```bash
module purge
module load intel_comp/19.1 impi/5.0.3 mumps/4.10.0 hdf5/1.8.9
```

1.3 Enter the `gyacomo` directory:
```bash
cd gyacomo
```

1.4 Make a `lib` directory, e.g., in `/gyacomo`:
```bash
mkdir lib
```

**Important:** Adjust the `gyacomo/local/dirs.inc` file by setting the paths to your `gyacomo` folder in the `PREFIX` variable and your path to the `lib` folder in the `LIBDIR` variable.

#### Test Installation and Linkage
Test the installation and linkage by typing `make` in `/gyacomo`. The compilation should halt with an error at the call of `futils` routines in `src/parallel_mod.F90`.

## Install Libraries
We now install the required external librairies.

### 2. Futils
2.1 Clone the `futils` library (or request access if needed):
```bash
git clone https://c4science.ch/diffusion/FUTILS/futils.git
```

2.2 Navigate to the `src` directory:
```bash
cd /futils/src
```

2.3 Adapt the makefile:
```bash
sed -i '37s/.*/OPT = -O3/' Makefile
```

2.4 Then compile:
```bash
make lib
```

2.5 Create necessary directories:
```bash
mkdir -p include/O
mkdir -p include/g
mkdir -p lib/O
mkdir -p lib/g
```

2.6 Install:
```bash
make install
```

#### Test Installation and Linkage
Test the installation and linkage by typing `make` in `/gyacomo`. The compilation should halt with an error at the call of `fftw` routines in `src.fourier_mod.F90`.

### 3. FFTW
3.1 Download `fftw` zip directory, unzip, and navigate to the directory:
```bash
wget http://www.fftw.org/fftw-3.3.10.tar.gz
tar -xvf fftw-3.3.10.tar.gz
cd fftw-3.3.10
```

3.2 Configure and make the double-precision version:
```bash
./configure --enable-mpi --prefix=$PWD
make
```

3.3 Install it:
```bash
make install
```

3.4 Repeat for the single-precision version:
```bash
./configure --enable-float --enable-mpi --prefix=$PWD
make
make install
```

3.5 Gather lib files in a new `lib` subdirectory:
```bash
mkdir lib
cp lib64/* lib/.
cp *.la lib/.
```

#### Navigate Back to Gyacomo
Go back to `/gyacomo` and verify that `FFTW3DIR   = $(LIBDIR)/fftw-3.3.10` in `/gyacomo/local/dirs.inc` when `ENVTYPE = Linux`.

#### Test Installation and Linkage
Test the installation and linkage by typing `make` in `/gyacomo`. The compilation should halt with an error at the call of `FM` routines in the `src/coeff_mod.F90` file.

## 4. FM
4.1 Navigate back to the main `lib` directory (e.g., `/gyacomo/lib/`) and download the zipped folder:
```bash
wget https://dmsmith.lmu.build/FM1.4/FM_files.zip
```

4.2 Unzip and enter the folder:
```bash
unzip FM_files.zip
mv FM_files FM
cd FM
```

4.3 Copy the source file with a `.F90` extension (for older Intel compilers):
```bash
for file in *.f95 ; do cp "$file" "${file%.*}.F90" ; done
```

4.4 Compile manually (~5min):
```bash
ifort fmsave.f95 -c -O3
ifort fm.f95 -c -O3
ifort fmzm90.f95 -c -O3
ifort TestFM.f95 -c -O3
ifort SampleFM.f95 -c -O3
ifort fmsave.o fm.o fmzm90.o TestFM.o -o TestFM
ifort fmsave.o fm.o fmzm90.o SampleFM.o -o SampleFM
```

4.5 Test the installation:
```bash
./TestFM
./SampleFM
```

4.6 Put the library together and move it to a local `lib` directory:
```bash
mkdir lib
ar r libfm.a fm.o fmsave.o fmzm90.o
mv libfm.a lib
```

4.7 Move the `.mod` files to a `mod` local directory:
```bash
mkdir mod
mv *.mod mod
```

## Test Installation and Linkage
Test the installation and linkage by typing `make` in `/gyacomo`. The compilation should proceed and produce an executable in `/gyacomo/bin/.` named `gyacomo23_dp`.

## 5. First Run
If the executable `gyacomo23_dp` is present in `/gyacomo/bin/.`, you can now test the first run of the code.
In `/gyacomo`, run the simulations setup script:
```bash
sh new_prob.sh ZBC
```
`ZBC`stands for Z-pinch base case, which is a 2D nonlinear simulation (recommended for a light run). You can also call the script with `CBC`, which is a minimal example of a cyclone base case (underresolved to be light).

This creates a new folder, `/gyacomo/simulations/`, with an example of a simulation directory `/gyacomo/simulations/problem_01`. Inside, you can find a `tutorial.md` file explaining how to run the code and analyze its data.