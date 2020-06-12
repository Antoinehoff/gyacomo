Hermite-Laguerre Z-pinch solver (2020)

Goal : Write the MOLI solver in Fortran and include a non linear Poisson bracket term to it.

Roadmap :
0) Write MOLI matlab solver in Fortran using Monli1D as starting point
	0.1) go from 1D space to 2D fourier and from Hermite basis to Hermite-Laguerre basis <- current stage
	0.2) implement linear Poisson equation in fourier space
	0.3) implement moment hierarchy linear terms
	0.4) implement explicit RK4 time solver
	0.5) Benchmark with MOLI matlab results for Z-pinch

1) Implementation of the non linear Poisson brackets term  
