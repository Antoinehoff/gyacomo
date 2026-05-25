"""
Example: Setup and Run a Gyacomo Simulation

This script demonstrates how to:
1. Create a new simulation from scratch
2. Modify parameters
3. Execute the simulation
4. Check status
"""

import sys
sys.path.append('..')
from gyacomo_simulation import GyacomoSimulation

# Create a new simulation instance
sim = GyacomoSimulation(
    run_dir='../../../simulations/example_run',
    jobnum=0
)

# Modify parameters
sim.params['BASIC']['tmax'] = 10.0
sim.params['BASIC']['dt'] = 0.001

# Grid parameters
sim.params['GRID']['Nx'] = 4
sim.params['GRID']['Ny'] = 32
sim.params['GRID']['Nz'] = 32
sim.params['GRID']['Lx'] = 200.0
sim.params['GRID']['Ly'] = 100.0

# Geometry (ITG-like case)
sim.params['GEOMETRY']['q0'] = 1.5
sim.params['GEOMETRY']['shear'] = 1.0
sim.params['GEOMETRY']['eps'] = 0.18
sim.params['GEOMETRY']['delta'] = 0.0  # No triangularity

# Model parameters
sim.params['MODEL']['LINEARITY'] = 'linear'
sim.params['MODEL']['nu'] = 0.05

# Species parameters (ITG unstable)
sim.params['SPECIES'][0]['k_N_'] = 2.0   # Ion density gradient
sim.params['SPECIES'][0]['k_T_'] = 6.0   # Ion temperature gradient
sim.params['SPECIES'][1]['k_N_'] = 2.0   # Electron density gradient
sim.params['SPECIES'][1]['k_T_'] = 2.0   # Electron temperature gradient

# Diagnostics
sim.params['DIAGNOSTICS']['dtsave_0d'] = 0.1
sim.params['DIAGNOSTICS']['dtsave_3d'] = 0.5

# Setup directory and write params
print("Setting up simulation directory...")
sim.setup_directory(clean=True)

print("\nSimulation parameters:")
print(f"  Run directory: {sim.run_dir}")
print(f"  Grid: Nx={sim.params['GRID']['Nx']}, Ny={sim.params['GRID']['Ny']}, Nz={sim.params['GRID']['Nz']}")
print(f"  Time: dt={sim.params['BASIC']['dt']}, tmax={sim.params['BASIC']['tmax']}")
print(f"  Executable: {sim.executable}")

# Uncomment to run the simulation
# print("\nRunning simulation...")
# result = sim.run(nproc=8, blocking=True, verbose=True)

# Check status
print("\nChecking status...")
status = sim.check_status()
print(f"  Output exists: {status['output_exists']}")
print(f"  Current time: {status['current_time']}")
print(f"  Progress: {status['progress']:.1f}%" if status['progress'] else "  Not started")
