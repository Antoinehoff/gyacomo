# HeLaZ (Hermite-Laguerre Z-pinch solver, 2021)
To compile it check INSTALLATION.txt

How to run it

1. Be sure to have correct paths in local/dirs.inc for the different libraries
2. Compile from HeLaZ/ using make
3. To run the code, use HeLaZ/wk/local_run.m and run it to set up the parameters and the results folder
4. Then go to the results folder and launch HeLaZ using mpirun -np num_procs ./../../../bin/helaz num_p num_kr
5. You can obtain various plots and gifs using HeLaZ/wk/analysis_2D.m once the simulation is done. To select the correct output file, run parameters*.m with the corresponding simulation parameters and then run analysis_2D.m (everything with matlab from wk/)

// Comment : For some collision operators (Sugama and Full Coulomb) you have to run COSOlver from B.J.Frei first in order to generate the required matrices in HeLaZ/iCa folder.

# Road map
(Current version : 2.5)

0. Write MOLI matlab solver in Fortran using Monli1D as starting point

	0.0 go from 1D space to 2D fourier and from Hermite basis to Hermite-Laguerre basis

	0.1 implement linear Poisson equation in fourier space

	0.2 implement moment hierarchy linear terms

	0.3 RK4 time solver

	0.4 Benchmark with MOLI matlab results for Z-pinch (cf. kz_linear script)

	0.5 Load COSOlver matrices

	0.6 Benchmarks now include Dougherty, Lenard-Bernstein and Full Coulomb collision operators

1. Implementation of the non linear Poisson brackets term

	1.0 FFTW3 has been used to treat the convolution as a product and discrete fourier transform

	1.1 Methods in fourier_mod.f90 have been validated by tests on Hasegawa Wakatani system

	1.1 Qualitative test : find similar turbulences as Hasegawa Wakatani system with few moments

	1.2 Zonal flows are observed in a similar way to Ricci Rogers 2006 with GS2

	1.3 Linear analysis showed that a certain amount of PJ are recquired to trigger mode

	1.4 Quantitative study with stationary average particle flux \Gamma_\infty

2. MPI parallel version

	2.1 First compilable parallel version (1D parallel along kr)

	2.2 Allow restart with different P,J values (results are not concluents)

	2.3 GK Dougherty operator

	2.4 2D cartesian parallel (along p and kr)

>2.5 GK Sugama collision operator

	2.6 GPU accelerated version

	2.7 GK Full Coulomb collision operator

4. GK 3D version, kr,kz,kpar for linear device

5. DK 3D version, kr,kz,kpar for linear device

6. DK+GK 3D version, kr,kz,kpar for linear device

7. 3D version with curvature
