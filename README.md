<figure>
<img src="https://c4scdn.ch/file/data/7a6vpqgtfcxtwhkpd4hu/PHID-FILE-wlsgn3omnbfilbqnzsvb/ezgif-2-ebfac79eeb26.gif" width="240">
<figcaption align = "center">
<i>Turbulence in a Z-pinch (the axis (r,z) corresponds to (x,y) in the current version of the code)</i>
</figcaption>
</figure>

GYACOMO (Gyrokinetic Advanced Collision Moment solver)
Copyright (C) 2022 EPFL

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

Author: Antoine C.D. Hoffmann

Contact: antoine.hoffmann@epfl.ch

##### Citing GYACOMO
If you use GYACOMO in your work, please cite (at least) one of the following paper: 

Hoffmann, A.C.D., Frei, B.J. & Ricci, P. (2023). Gyrokinetic simulations of plasma turbulence in a Z-pinch using a moment-based approach and advanced collision operators. Journal of Plasma Physics, 89(2), 905890214. doi:10.1017/S0022377823000284

Hoffmann, A.C.D., Frei, B.J. & Ricci, P. (2023). Gyrokinetic moment-based simulations of the Dimits shift. Journal of Plasma Physics, 89(6), 905890611. doi:10.1017/S0022377823001320

# What is GYACOMO ?

GYACOMO is the Gyrokinetic Advanced Collision Moment solver which solves the gyrokinetic Boltzmann equation in the delta-f flux-tube limit based on a projection of the velocity distribution function onto a Hermite-Laguerre velocity basis.

It can be coupled with precomputed matrices from the code Cosolver (B.J. Frei) to incorporate advanced collision operators up to the gyro-averaged linearized exact coulomb interaction (GK Landau operator).

This repository contains the solver source code (in /src) but also my personnal post-processing Matlab scripts, which are less documented. I would recommend the user to write their own post-processing scripts based on the H5 files the code outputs.

#### GYACOMO can
- run in parallel using MPI (mpirun -np N ./path_to_exec Np Ny Nz, where N = Np x Ny x Nz is the number of processes and Np Ny Nz are the parallel dimensions in Hermite polynomials, binormal direction, and parallel direction, respectively).
- run in single precision.
- evolve kinetic electrons and ions.
- use an adiabatic electrons model.
- include perpendicular magnetic fluctuations.
- use Z-pinch, s-alpha, circular and Miller geometry model.
- use various experimental closures for the linear and nonlinear terms.
- use linear GK Landau, Sugama, Lorentz collision operators. (requires precomputed matrix files, ask them!)
- add background ExB shear flow. (Hammett's method)
- use an adiabatic ion model. (not verified)
#### GYACOMO cannot (I wish it could...)
- include parallel magnetic field fluctuations. (easy)
- include finite rhostar effects. (hard)
- run without the futils library. (easy but boring, ask the zip file!)
- Use shared memory parallelization. (okish)
- run global simulations. (for another code)

# How to compile and run GYACOMO

A tutorial is now on the code's wiki https://gitlab.epfl.ch/ahoffman/gyacomo/-/wikis/home.

(shorter guideline)
To compile and run GYACOMO, follow these steps:

1. Make sure the correct library paths are set in local/dirs.inc. Refer to INSTALATION.txt for instructions on installing the required libraries.
2. Go to the /gyacomo directory and run make to compile the code. The resulting binary will be located in /gyacomo/bin. You can also compile a debug version by running make dbg.
3. Once the compilation is done, you can test the executable in the directory /testcases/cyclone_example (see README there) or /testcases/zpinch_example (see wiki).
4. To stop the simulation without corrupting the output file, create a blank file called "mystop" using $touch mystop in the directory where the simulation is running. The code looks for this file every 100 time steps and the file will be removed once the simulation is stopped.
5. It is possible to chain simulations by using the parameter "Job2load" in the fort_XX.90 file. For example, to restart a simulation from the latest 5D state saved in outputs_00.h5, create a new fort.90 file called fort_01.90 and set "Job2load" to 0. Then run GYACOMO with the command ./gyacomo 0 or mpirun -np N ./gyacomo Np Ny Nz 0. This will create a new output file called output_01.h5.
6. To generate plots and gifs using the simulation results, use the script gyacomo/wk/gyacomo_analysis.m and specify the directory where the results are located. Note that this script is not currently a function.

Note: For some collision operators (Sugama and Full Coulomb), you will need to run COSOlver from B.J.Frei to generate the required matrices in the gyacomo/iCa folder before running GYACOMO.



# Changelog

### 4.x GYACOMO

>4.11 background ExB shear

>4.1 Miller geometry benchmarked

>4.01 Singular value decomposition is availale with LAPACK (used for DLRA experiments)

>4.0 Gyacomo is born and the code is open-source with a GNU GPLv3 license

### 3.x HeLaZ 3D (flux tube s-alpha)

>3.9 Perpendicular electromagnetic fluctuations by solving Ampere equations (benchmarked linearly)

>3.8 Benchmarked for CBC against GENE for various gradients values (see Dimits_fig3.m)

>3.7 The frequency plane has been transposed from positive kx to positive ky for easier implementation of shear. Also added 3D Z-pinch geometry

>3.6 MPI 3D parallelization in p, kx and z and benchmarked for each parallel options with gbms (new molix) for linear fluxtube shearless.

>3.5 Staggered grid for parallel odd/even coupling

>3.4 Adiabatic electrons

>3.3 Benchmarked in fluxtube s-alpha geometry linear run with molix (B.J.Frei) code and works now for shear = 0 with periodic z BC

>3.2 Stopping file procedure like in GBS is added

>3.1 Implementation of mirror force

>3.0 3D version and works as the 2D version if Nz = 1, the coordinates were renamed from (r,z)  to (x,y,z). Now the parallel direction is ez.

### 2.x 2D Zpinch MPI parallel version

>2.7 Versatile interpolation of kperp for the cosolver matrices and corrections done on DGGK

>2.6 Change of collisionality normalisation (from nu_ei to nu_ii), implementation of FCGK

>2.5 GK cosolver collision implementation

>2.4 MPI 2D cartesian parallel (along p and kr)

>2.3 GK Dougherty operator

>2.2 Allow restart with different P,J values

>2.1 First compilable parallel version (1D parallel along kr)

### 1.x Implementation of the non linear Poisson bracket term

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
