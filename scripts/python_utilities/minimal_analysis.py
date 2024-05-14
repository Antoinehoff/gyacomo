"""
This script loads data from the HDF5 output file of GYACOMO ('outputs_00.h5').
It plots the time traces of the heat and particle fluxes in the radial direction by loading
the following datasets:
- 'data/var0d/time': Time values
- 'data/var0d/hflux_x': Heat flux data
- 'data/var0d/pflux_x': Particle flux data
"""
import load_data as loader
import matplotlib.pyplot as plt
import numpy as np
import sys
import tools


#------ Open the HDF5 file
if len(sys.argv) == 1:
    # No command-line argument provided, use default filename
   filename = "outputs_00.h5"
elif len(sys.argv) == 2:
    # Command-line argument provided
    arg = sys.argv[1]
    if arg.isdigit():
        # Argument is an integer, format it as XX
        filename = f"outputs_{int(arg):02d}.h5"
    else:
        # Argument is a string
        filename = arg
else:
    print("Usage: python minimal_analysis.py [filename, int or string (opt)]")
    sys.exit(1)

#------ Load data
# load the input parameters from the STDIN.00 file stored in the h5
params      = loader.load_params(filename)

# time traces (heat and particle fluxes)
t0D, hflux_x = loader.load_data_0D(filename,'hflux_x')
t0D, pflux_x = loader.load_data_0D(filename,'pflux_x')
# 3D fields (phi and first gyromoment)
t3D, phi, tf = loader.load_data_3D_frame(filename,'phi',10000)
t3D, Ni00,tf = loader.load_data_3D_frame(filename,'Na00',10000)
# grids of the simulation
x,kx,y,ky,z,p,j  = loader.load_grids(filename)

# Examples to load data in dictionaries
# metric (gij, dBdi, Jacobian etc.)
metric      = loader.load_group(filename,'metric')
# another way to load the time traces all at once
time_traces = loader.load_group(filename,'var0d')

#------ Process data
# Fourier transform into real space at constant plane
phi = tools.zkxky_to_xy_const_z(phi,-1)
Ni00= tools.zkxky_to_xy_const_z(Ni00,-1)

#------ Plot data
# Plot hflux_x and pflux_x against time
fig, axes = plt.subplots(2, 2, figsize=(10, 6))
# Plot hflux_x
nt0D = t0D.size
if hflux_x.size > nt0D :
    axes[0,0].plot(t0D, hflux_x[::2], label='ion heat flux')
    axes[0,0].plot(t0D, hflux_x[1::2], label='electron heat flux')
else:
    axes[0,0].plot(t0D, hflux_x, label='ion heat flux')
axes[0,0].set_title('radial heat flux')
axes[0,0].set_xlabel('t c_s/R')
axes[0,0].set_ylabel('Q_x')

# Plot pflux_x
if pflux_x.size > nt0D :
    axes[1,0].plot(t0D, pflux_x[::2], label='ion heat flux')
    axes[1,0].plot(t0D, pflux_x[1::2], label='electron heat flux')
else:
    axes[0,0].plot(t0D, pflux_x, label='ion heat flux')
axes[1,0].set_title('radial particle flux')
axes[1,0].set_xlabel('t c_s/R')
axes[1,0].set_ylabel('P_x')

# Plot slice of electrostatic potential at constant z
axes[0,1].imshow(Ni00, extent=[x[0], x[-1], y[0], y[-1]],cmap='viridis',interpolation='quadric')
#axes[0,1].pcolor(x,y,Ni00, cmap='viridis')
axes[0,1].set_title(f'Ni00 (z = 0, t = {tf:3.1f})')
axes[0,1].set_xlabel('x')
axes[0,1].set_ylabel('y')
#fig.colorbar(im, ax=axes[0,1], label='Ni00')

# Plot slice of electrostatic potential at constant z
axes[1,1].imshow(phi, extent=[x[0], x[-1], y[0], y[-1]],cmap='viridis',interpolation='quadric')
#axes[1,1].pcolor(x,y,phi, cmap='viridis')
axes[1,1].set_title(f'phi (z = 0, t = {tf:3.1f})')
axes[1,1].set_xlabel('x')
axes[1,1].set_ylabel('y')
#fig.colorbar(im, ax=axes[1,1], label='Potential')

plt.tight_layout()
plt.show()
