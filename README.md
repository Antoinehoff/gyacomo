# HeLaZ (Hermite-Laguerre Z-pinch solver, 2021)
To compile it check INSTALLATION.txt

How to run it

1. Be sure to have correct paths in local/dirs.inc for the different libraries
2. Compile from HeLaZ/ using make
3. To run the code, use HeLaZ/wk/local_run.m and run it to set up the parameters and the results folder
4. Then go to the results folder and launch HeLaZ using mpirun -np num_procs ./../../../bin/helaz num_p num_kr
5. You can obtain various plots and gifs using HeLaZ/wk/analysis_2D.m once the simulation is done. To select the correct output file, run parameters*.m with the corresponding simulation parameters and then run analysis_2D.m (everything with matlab from wk/)
6. The current simulation can be stopped at any moment by writing a file named "stop" in the simulation directory
// Comment : For some collision operators (Sugama and Full Coulomb) you have to run COSOlver from B.J.Frei first in order to generate the required matrices in HeLaZ/iCa folder.

# Changelog

3. HeLaZ 3D

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

# Roadmap

3. 3D version

	3.1 Flux-tube

	3.2 Benchmark and 3D transport analysis

4. full-f terms? Evolving backgrounds?
