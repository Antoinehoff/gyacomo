"""
Example: Parameter Scan

This script demonstrates how to:
1. Set up multiple simulations with varying parameters
2. Run parameter scans (e.g., density gradient scan)
3. Organize results in separate directories
"""

import sys
sys.path.append('..')
from gyacomo_simulation import GyacomoSimulation
import numpy as np

# Define parameter scan
base_dir = '../../../simulations/parameter_scan'
k_N_values = [1.0, 2.0, 4.0, 6.0, 8.0]  # Density gradient scan

print("Setting up parameter scan...")
print(f"Base directory: {base_dir}")
print(f"Scanning k_N from {k_N_values[0]} to {k_N_values[-1]}")
print(f"Number of cases: {len(k_N_values)}")

simulations = []

for i, k_N in enumerate(k_N_values):
    # Create simulation in subdirectory
    run_dir = f"{base_dir}/k_N_{k_N:.1f}"
    
    print(f"\nCase {i+1}: k_N = {k_N:.1f}")
    print(f"  Directory: {run_dir}")
    
    # Create simulation instance
    sim = GyacomoSimulation(run_dir=run_dir, jobnum=0)
    
    # Set base parameters
    sim.params['BASIC']['tmax'] = 10.0
    sim.params['BASIC']['dt'] = 0.001
    
    # Grid
    sim.params['GRID']['Nx'] = 2
    sim.params['GRID']['Ny'] = 24
    sim.params['GRID']['Nz'] = 24
    
    # Scan parameter: density gradient
    sim.params['SPECIES'][0]['k_N_'] = k_N  # Ion
    sim.params['SPECIES'][1]['k_N_'] = k_N  # Electron
    
    # Fixed temperature gradient (ITG-like)
    sim.params['SPECIES'][0]['k_T_'] = 6.0
    sim.params['SPECIES'][1]['k_T_'] = 2.0
    
    # Setup directory
    sim.setup_directory(clean=True)
    print(f"  Created directory and wrote params.in")
    
    simulations.append(sim)

print("\n" + "="*60)
print("SETUP COMPLETE")
print("="*60)
print(f"Created {len(simulations)} simulation directories")
print("\nTo run all simulations, you can:")
print("  1. Manually cd to each directory and run:")
print("     mpirun -np 8 /path/to/gyacomo.exe 1 8 1 < params.in")
print("  2. Use a batch submission script")
print("  3. Uncomment the following loop to run sequentially:\n")

# Uncomment to run all simulations sequentially
# print("Running simulations...")
# for i, sim in enumerate(simulations):
#     print(f"\nRunning case {i+1}/{len(simulations)}...")
#     result = sim.run(nproc=8, blocking=True, verbose=True)
#     status = sim.check_status()
#     print(f"  Completed: {status['progress']:.1f}%")

print("\nExample analysis after runs complete:")
print("""
import matplotlib.pyplot as plt
import numpy as np

# Collect results
growth_rates = []
for sim in simulations:
    result = sim.compute_growth_rate(field_name='phi', kx_val=0.0)
    gamma_max = np.real(result['omega']).max()
    growth_rates.append(gamma_max)

# Plot scan results
plt.figure(figsize=(8, 6))
plt.plot(k_N_values, growth_rates, 'o-')
plt.xlabel(r'$k_N$ [density gradient]')
plt.ylabel(r'$\gamma_{max}$ [c_s/R]')
plt.title('Maximum Growth Rate vs Density Gradient')
plt.grid(True, alpha=0.3)
plt.savefig('k_N_scan.png', dpi=150)
plt.show()
""")
