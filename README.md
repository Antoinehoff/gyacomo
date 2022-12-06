![Alt Text](https://c4scdn.ch/file/data/7a6vpqgtfcxtwhkpd4hu/PHID-FILE-wlsgn3omnbfilbqnzsvb/ezgif-2-ebfac79eeb26.gif)
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

1. Be sure to have correct library paths in local/dirs.inc for the different libraries, see INSTALATION.txt for a tutorial to install the required libraries.
2. Compile from /gyacomo using make, the binary will be located in /gyacomo/bin (you can also compile a debug version using make dbg)
4. The fort.90 file should contain the parameters for a typical CBC to test the compilation. One can run it by calling the executable /bin/gyacomo in the directory where the fort.90 is located.
5. It is possible to run GYACOMO in parallel using MPI: mpirun -np N ./bin/gyacomo Np Ny Nz where N=Np x Ny x Nz is the number of processes and Np Ny Nz are the parallel dimensions in Hermite polynomials, binormal direction and parallel direction, respectively
6. You can stop your simulation without corrupting the output file by creating a blank file call "mystop", using e.g. "touch mystop" in the directory where the simulation is running. (the file will be removed once read)
7. It is also possible to put simulations ID in order to chain them. The parameter Job2load allows you to tell which output file should be read in order to restart a simulation. E.g. I run a first simulation with job2load = -1, it creates a outputs_00.h5 then I create a new fort.90 which I call fort_01.90 where job2load = 0. I run then GYACOMO, indicating that I want it to read the fort_00.90 using 0 as a last argument, i.e. "./gyacomo 0" or "mpirun -np N ./gyacomo Np Ny Nz 0", which will start from the latest 5D state saved in outputs_00.h5. A new output file has also been created, output_01.h5.
8. You can obtain various plots and gifs using gyacomo/wk/gyacomo_analysis.m once the simulation is done. The directory where the results are located must be given in the scripts. It is not a function (yet...)
// Comment : For some collision operators (Sugama and Full Coulomb) you have to run COSOlver from B.J.Frei first in order to generate the required matrices in gyacomo/iCa folder. //



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
