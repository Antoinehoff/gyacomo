# Gyacomo (Gyrokinetic Advanced Collision Moment Solver)

![ezgif-2-ebfac79eeb26](https://github.com/user-attachments/assets/e38bbeed-e672-4a32-a6c9-2e321086656b)

*Fig. 1: Turbulence and zonal flows in a Z-pinch.*

---

**Gyacomo** is a solver for the gyrokinetic Boltzmann equation in the delta-f flux-tube limit. It is based on a Hermite-Laguerre moment expansion of the velocity distribution function. The code supports kinetic and adiabatic electrons/ions, multiple geometries, advanced collision operators, and parallelization via MPI.

Originally developed at EPFL, Gyacomo is now open-source under the GNU General Public License v3.

**Author**: Antoine C.D. Hoffmann  
**Contact**: ahoffmann@pppl.gov  

---

## 📄 License

This program is free software: you can redistribute it and/or modify it under the terms of the [GNU General Public License](https://www.gnu.org/licenses/) as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

It is distributed in the hope that it will be useful, but **WITHOUT ANY WARRANTY**—without even the implied warranty of **MERCHANTABILITY** or **FITNESS FOR A PARTICULAR PURPOSE**.

---

## 📚 References

If you use or mention **Gyacomo** in your work, please cite the following reference:

- Hoffmann, A.C.D., Frei, B.J. & Ricci, P. (2023). *Gyrokinetic moment-based simulations of the Dimits shift.* J. Plasma Phys. 89(6), 905890611.  
  [doi:10.1017/S0022377823001320](https://doi.org/10.1017/S0022377823001320)

Publications:

- Hoffmann, A.C.D., Frei, B.J. & Ricci, P. (2023). *Gyrokinetic simulations of plasma turbulence in a Z-pinch using a moment-based approach and advanced collision operators.* J. Plasma Phys. 89(2), 905890214.  
  [doi:10.1017/S0022377823000284](https://doi.org/10.1017/S0022377823000284)

- Hoffmann, A.C.D., Frei, B.J. & Ricci, P. (2023). *Gyrokinetic moment-based simulations of the Dimits shift.* J. Plasma Phys. 89(6), 905890611.  
  [doi:10.1017/S0022377823001320](https://doi.org/10.1017/S0022377823001320)

- Frei, B.J. (2023). *A Gyrokinetic Moment Model of the Plasma Boundary in Fusion Devices.* EPFL, Lausanne.   
  [doi:10.5075/epfl-thesis-9960](https://doi.org/10.5075/epfl-thesis-9960)

- Hoffmann, A.C.D., Frei, B., Giroud-Garampon, P., & Ricci, P. (2024). *A gyrokinetic moment-based approach for multi-scale multi-fidelity turbulence simulations*. [Computational Challenges and Optimization in Kinetic Plasma Physics](https://www.imsi.institute/computational-challenges-and-optimization-in-kinetic-plasma-physics-poster-session/), Institute for Mathematical Science and Statistics Innovation, Chicago, IL.

- Hoffmann, A.C.D., Balestri, A., & Ricci, P. (2025). *Investigation of triangularity effects on tokamak edge turbulence through multi-fidelity gyrokinetic simulations.* Plasma Phys. Control. Fusion 67 015031.   
  [doi:10.1088/1361-6587/ad9e6f](https://iopscience.iop.org/article/10.1088/1361-6587/ad9e6f)

- Hoffmann, A.C.D. (2025). *Nonlinear Simulation of Plasma Turbulence Using a Gyrokinetic Moment-Based Approach.* EPFL, Lausanne, 2024.   
[doi:10.5075/epfl-thesis-10651](https://doi.org/10.5075/epfl-thesis-10651)

---

## ❓ What is Gyacomo?

Gyacomo solves the gyrokinetic Boltzmann equation projected onto a Hermite-Laguerre velocity space basis. It supports flux-tube simulations in multiple geometries and includes kinetic effects, electromagnetic fluctuations, and advanced collisions using precomputed matrices from **Cosolver** (by B.J. Frei).

> 💡 This repository also includes personal post-processing MATLAB scripts. They're provided as-is and less documented — users are encouraged to write their own post-processing routines using the HDF5 output files.

---

## ✅ Gyacomo Can

- Run in parallel using MPI:  
  `mpirun -np N ./path_to_exec Np Ny Nz`  
  where `N = Np × Ny × Nz` (parallelization in Hermite modes, binormal, and parallel directions).

- Run in **single precision**

- Handle **kinetic or adiabatic electrons and ions**

- Include **perpendicular magnetic fluctuations**

- Simulate in **Z-pinch**, **s-alpha**, **circular**, and **Miller** geometries

- Use various closures for linear and nonlinear terms

- Incorporate linear GK Landau, Sugama, and Lorentz collision operators  
  *(requires precomputed matrices from Cosolver — ask for them!)*

- Add background **ExB shear flows** (Hammett’s method)

- Use an **adiabatic ion model** *(not fully verified)*

---

## ❌ Gyacomo Cannot (yet)

- Include **parallel magnetic field fluctuations** *(should be easy)*

- Model **finite ρ\* (rhostar) effects** *(hard)*

- Run without the **futils** library *(easy, but tedious — ask for the zip!)*

- Use **shared memory** parallelization *(semi-working)*

- Run **global simulations** *(out of scope for now)*

---

## ⚙️ How to Compile and Run

Detailed installation and usage instructions are available on the [Gyacomo Wiki](https://gitlab.epfl.ch/ahoffman/Gyacomo/-/wikis/home).

> 🔧 Note: Some collision operators (Sugama, full Coulomb) require you to precompute matrix files using **COSOlver** and place them in the `Gyacomo/iCa` folder before running simulations.

---

## 🕒 Changelog

### v3.x – Gyacomo

- Initial public release with [GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)
- Python post-processing utilities
- Installation tutorials
- Background **ExB shear** implemented
- **Miller geometry** benchmarked
- Support for **SVD** with LAPACK (for DLRA experiments)

---

### v2.x – HeLaZ 3D (Flux-tube s-alpha)

- Perpendicular electromagnetic fluctuations (Ampère’s law, linearly benchmarked)
- Benchmarked for CBC against GENE (see `Dimits_fig3.m`)
- Sheared implementation via frequency transpose (kx → ky)
- Full 3D MPI parallelization (p, kx, z)
- **Staggered grid** for parallel odd/even coupling
- **Adiabatic electron** model
- Mirror force implementation
- 2D/3D hybrid mode (Nz=1 = 2D)
- Restart capability with different P, J
- GBS-style stop file procedure

---

### v1.x – 2D Z-pinch MPI Parallel Version

- Parallelization in 2D (p, kr)
- Implementation of GK Dougherty and Cosolver-based collision operators
- New collisionality normalization (from νₑᵢ to νᵢᵢ)
- Versatile k_perp interpolation for Cosolver matrices
- Compatibility with MOLI MATLAB benchmarks
- First stable parallel release

---

### v0.x – Initial Implementation

- Nonlinear Poisson bracket term
- RK4 time solver
- Moment-hierarchy linear terms
- Fourier space linear Poisson equation
- Validation against Hasegawa-Wakatani system
- FFTW3 for convolution
- Support for Hermite and Hermite-Laguerre bases
- Start from GBS skeleton

---

Feel free to reach out for matrix files, help compiling, or contributing to the project!
