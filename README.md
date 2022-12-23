<figure>
<img src="https://c4scdn.ch/file/data/7a6vpqgtfcxtwhkpd4hu/PHID-FILE-wlsgn3omnbfilbqnzsvb/ezgif-2-ebfac79eeb26.gif" width="240">
<figcaption align = "center">
<i>Turbulence in a Z-pinch (the axis (r,z) corresponds to (x,y) in the current version of the code)</i>
</figcaption>
</figure>

GYACOMO (Gyrokinetic Advanced Collision Moment solver, 2021)
Copyright (C) 2022  A.C.D. Hoffmann

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

# How to compile and run GYACOMO

To compile and run GYACOMO, follow these steps:


1. Make sure the correct library paths are set in local/dirs.inc. Refer to INSTALATION.txt for instructions on installing the required libraries.
2. Go to the /gyacomo directory and run make to compile the code. The resulting binary will be located in /gyacomo/bin. You can also compile a debug version by running make dbg.
3. The file fort.90 should contain the parameters for a typical CBC. To test the compilation, navigate to the directory where fort.90 is located and run the executable /bin/gyacomo.
4. GYACOMO can be run in parallel using MPI by running mpirun -np N ./bin/gyacomo Np Ny Nz, where N = Np x Ny x Nz is the number of processes and Np Ny Nz are the parallel dimensions in Hermite polynomials, binormal direction, and parallel direction, respectively.
5. To stop the simulation without corrupting the output file, create a blank file called "mystop" using touch mystop in the directory where the simulation is running. The file will be removed once it is read.
6. It is possible to chain simulations by using the parameter "Job2load" in the fort.90 file. For example, to restart a simulation from the latest 5D state saved in outputs_00.h5, create a new fort.90 file called fort_01.90 and set "Job2load" to 0. Then run GYACOMO with the command ./gyacomo 0 or mpirun -np N ./gyacomo Np Ny Nz 0. This will create a new output file called output_01.h5.
7. To generate plots and gifs using the simulation results, use the script gyacomo/wk/gyacomo_analysis.m and specify the directory where the results are located. Note that this script is not currently a function.
Note: For some collision operators (Sugama and Full Coulomb), you will need to run COSOlver from B.J.Frei to generate the required matrices in the gyacomo/iCa folder before running GYACOMO.



# Changelog

### 4.x GYACOMO

>4.1 Miller geometry is added and benchmarked for CBC adiabatic electrons

>4.0 new name and opening the code with GNU GPLv3 license

### 3.x HeLaZ 3D (flux tube s-alpha)

>3.9 HeLaZ can now evolve electromagnetic fluctuations by solving Ampere equations (benchmarked linearly)

>3.8 HeLaZ has been benchmarked for CBC with GENE for various gradients values (see Dimits_fig3.m)

>3.7 The frequency plane has been transposed from positive kx to positive ky for easier implementation of shear. Also added 3D zpinch geometry

>3.6 HeLaZ is now parallelized in p, kx and z and benchmarked for each parallel options with gbms (new molix) for linear fluxtube shearless.

>3.5 Staggered grid for parallel odd/even coupling

>3.4 HeLaZ can run with adiabatic electrons now!

>3.3 HeLaZ 3D has been benchmarked in fluxtube salphaB geometry linear run with molix (B.J.Frei) code and works now for shear = 0 with periodic z BC

>3.2 Stopping file procedure like in GBS is added

>3.1 Implementation of mirror force

>3.0 HeLaZ is now 3D and works like HeLaZ 2D if Nz = 1, the axis were renamed from r and z  to x,y and z. Now the parallel direction is ez. All arrays have been extended, diagnostics and analysis too. The linear coefficients are now precomputed with geometry routines.

### 2.x 2D Zpinch MPI parallel version

>2.7 Versatile interpolation of kperp for the cosolver matrices and corrections done on DGGK

>2.6 Change of collisionality normalisation (from nu_ei to nu_ii), implementation of FCGK

>2.5 GK cosolver collision implementation

>2.4 2D cartesian parallel (along p and kr)

>2.3 GK Dougherty operator

>2.2 Allow restart with different P,J values (results are not concluents)

>2.1 First compilable parallel version (1D parallel along kr)

### 1.x Implementation of the non linear Poisson brackets term

>1.4 Quantitative study with stationary average particle flux \Gamma_\infty

>1.3 Linear analysis showed that a certain amount of PJ are recquired to trigger mode

>1.2 Zonal flows are observed in a similar way to Ricci Rogers 2006 with GS2

>1.1 Qualitative test : find similar turbulences as Hasegawa Wakatani system with few moments

>1.1 Methods in fourier_mod.f90 have been validated by tests on Hasegawa Wakatani system

>1.1 Methods in fourier_mod.f90 have been validated by tests on Hasegawa Wakatani system

>1.0 FFTW3 has been used to treat the convolution as a product and discrete fourier transform

### 0.x Write MOLI matlab solver in Fortran using Monli1D as starting point

>0.6 Benchmarks now include Dougherty, Lenard-Bernstein and Full Coulomb collision operators

>0.5 Load COSOlver matrices

>0.4 Benchmark with MOLI matlab results for Z-pinch (cf. kz_linear script)

>0.3 RK4 time solver

>0.2 implement moment hierarchy linear terms

>0.1 implement linear Poisson equation in fourier space

>0.0 go from 1D space to 2D fourier and from Hermite basis to Hermite-Laguerre basis
