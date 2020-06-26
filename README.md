HeLaZ (Hermite-Laguerre Z-pinch solver, 2020)

Roadmap :

0. Write MOLI matlab solver in Fortran using Monli1D as starting point

	0.1 go from 1D space to 2D fourier and from Hermite basis to Hermite-Laguerre basis

	0.2 implement linear Poisson equation in fourier space

	0.3 implement moment hierarchy linear terms

	0.4 RK4 time solver

	0.5 Benchmark with MOLI matlab results for Z-pinch (cf. kz_linear script)
		Note : benchmark.m points out differences between MOLI and HeLaZ 
		       in extremely simple cases as P=J=2, eta_X = 0 except one.

	0.6 Load COSOlver matrices <- Current stage

1. Implementation of the non linear Poisson brackets term  

