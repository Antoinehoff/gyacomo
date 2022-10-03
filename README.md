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
To compile it check INSTALLATION.txt

How to run it

1. Be sure to have correct library paths in local/dirs.inc for the different libraries
2. Compile from /gyacomo using make, the binary will be located in /gyacomo/bin
4. You can run a typical CBC to test the compilation using the basic fort.90 parameter file,
   just type ./bin/gyacomo
5. It is possible to run it in parallel (MPI) as mpirun -np N ./bin/gyacomo Np Ny Nz
   where N=Np x Ny x Nz is the number of processes and Np Ny Nz are the parallel dimensions in
	 Hermite polynomials, binormal direction and parallel direction, respectively
6. You can stop your simulation without breaking output file by creating a blank file call "mystop"
   in the directory where the simulation is running. (the file will be removed once read)
7. You can obtain various plots and gifs using gyacomo/wk/header_3D_results.m once the simulation is done.
// Comment : For some collision operators (Sugama and Full Coulomb) you have to run COSOlver from B.J.Frei first in order to generate the required matrices in gyacomo/iCa folder.

# Changelog
4. GYACOMO
  4.0 new naming and opening the code with GNU GPLv3 license

3. HeLaZ 3D
  3.9 HeLaZ can now evolve electromagnetic fluctuations by solving Ampere equations (benchmarked linearly)

	3.8 HeLaZ has been benchmarked for CBC with GENE for various gradients values (see Dimits_fig3.m)

	3.7 The frequency plane has been transposed from positive kx to positive ky for easier implementation of shear. Also added 3D zpinch geometry

	3.6 HeLaZ is now parallelized in p, kx and z and benchmarked for each parallel options with gbms (new molix) for linear fluxtube shearless.

	3.5 Staggered grid for parallel odd/even coupling

	3.4 HeLaZ can run with adiabatic electrons now!

	3.3 HeLaZ 3D has been benchmarked in fluxtube salphaB geometry linear run with molix (B.J.Frei) code and works now for shear = 0 with periodic z BC

	3.2 Stopping file procedure like in GBS is added

	3.1 Implementation of mirror force

	3.0 HeLaZ is now 3D and works like HeLaZ 2D if Nz = 1, the axis were renamed (r,z) -> (x,y,z) and now the parallel direction is ez. All arrays have been extended, diagnostics and analysis too. The linear coefficients are now precomputed with lin_coeff_and_geometry routines.

2. MPI parallel version

	2.7 Versatile interpolation of kperp for the cosolver matrices and corrections done on DGGK

	2.6 Change of collisionality normalisation (from nu_ei to nu_ii), implementation of FCGK

	2.5 GK cosolver collision implementation

	2.4 2D cartesian parallel (along p and kr)

	2.3 GK Dougherty operator

	2.2 Allow restart with different P,J values (results are not concluents)

	2.1 First compilable parallel version (1D parallel along kr)

1. Implementation of the non linear Poisson brackets term

	1.4 Quantitative study with stationary average particle flux \Gamma_\infty

	1.3 Linear analysis showed that a certain amount of PJ are recquired to trigger mode

	1.2 Zonal flows are observed in a similar way to Ricci Rogers 2006 with GS2

	1.1 Qualitative test : find similar turbulences as Hasegawa Wakatani system with few moments

	1.1 Methods in fourier_mod.f90 have been validated by tests on Hasegawa Wakatani system

	1.1 Methods in fourier_mod.f90 have been validated by tests on Hasegawa Wakatani system

	1.0 FFTW3 has been used to treat the convolution as a product and discrete fourier transform

0. Write MOLI matlab solver in Fortran using Monli1D as starting point

	0.6 Benchmarks now include Dougherty, Lenard-Bernstein and Full Coulomb collision operators

	0.5 Load COSOlver matrices

	0.4 Benchmark with MOLI matlab results for Z-pinch (cf. kz_linear script)

	0.3 RK4 time solver

	0.2 implement moment hierarchy linear terms

	0.1 implement linear Poisson equation in fourier space

	0.0 go from 1D space to 2D fourier and from Hermite basis to Hermite-Laguerre basis
