HeLaZ (Hermite-Laguerre Z-pinch solver, 2020)

Roadmap :

0. Write MOLI matlab solver in Fortran using Monli1D as starting point

	0.1 go from 1D space to 2D fourier and from Hermite basis to Hermite-Laguerre basis (Done)

	0.2 implement linear Poisson equation in fourier space (Done)

	0.3 implement moment hierarchy linear terms (Done)

	0.4 RK4 time solver (Done)

	0.5 Benchmark with MOLI matlab results for Z-pinch <- Current stage

1. Implementation of the non linear Poisson brackets term  

