# HeLaZ (Hermite-Laguerre Z-pinch solver, 2020)

How to run it :

1. Be sure to have correct paths in local/dirs.inc for the different libraries
2. You can compile from HeLaZ/ using make and launch from HeLaZ/wk using ./../bin/HeLaZ
3. To have a better interface, open a script HeLaZ/wk/parameters*.m and run it to set up a wanted simulation.
4. You can obtain various plots and gifs using HeLaZ/wk/analysis_2D.m once the simulation is done. To select the correct output file, run parameters*.m with the corresponding simulation parameters and then run analysis_2D.m (everything with matlab from wk/)

Roadmap : (Current version 1.4)

0. Write MOLI matlab solver in Fortran using Monli1D as starting point

	0.0 go from 1D space to 2D fourier and from Hermite basis to Hermite-Laguerre basis

	0.1 implement linear Poisson equation in fourier space

	0.2 implement moment hierarchy linear terms

	0.3 RK4 time solver

	0.4 Benchmark with MOLI matlab results for Z-pinch (cf. kz_linear script)
		Note : benchmark_*.m compares MOLI and HeLaZ linear results

	0.5 Load COSOlver matrices

	0.6 Benchmarks now include Dougherty, Lenard-Bernstein and Full Coulomb collision operators
	    Note : for full Coulomb, one must store a precomputed matrix from COSOlver in the iCa folder

1. Implementation of the non linear Poisson brackets term

	1.0 FFTW3 has been used to treat the convolution as a product and discrete fourier transform

	1.1 Methods in fourier_mod.f90 have been validated by tests on Hasegawa Wakatani system

	1.1 Qualitative test : find similar turbulences as Hasegawa Wakatani system with few moments

	1.2 Zonal flows are observed in a similar way to Ricci Rogers 2006 with GS2

	1.3 Linear analysis showed that a certain amount of PJ are recquired to trigger mode

		1.3.1 The \eta_B = 0.5 case is easier since it converged better in linear analysis than \eta_B = 1.0

		1.3.2 Collisionality helps

	1.4 Quantitative study with stationary average particle flux \Gamma_\infty

		1.4.1 Convergence study of \Gamma_\infty w.r.t. P and J

		1.4.2 Direct comparison with GS2 results of Ricci,Rogers 2006
