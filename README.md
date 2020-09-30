HeLaZ (Hermite-Laguerre Z-pinch solver, 2020)
Current version : 0.6

Roadmap :

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
