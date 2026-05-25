"""
Example: Analyze Gyacomo Simulation Results

This script demonstrates how to:
1. Load an existing simulation
2. Extract time traces
3. Load 3D field data
4. Create visualizations
"""

import sys
sys.path.append('..')
from gyacomo_simulation import GyacomoSimulation
import matplotlib.pyplot as plt
import numpy as np

# Load an existing simulation
# Replace with actual simulation directory
sim = GyacomoSimulation(
    run_dir='../../../simulations/check_instability_tcv_gkyl/PT_rho0.9',
    jobnum=0
)

# Check if output exists
status = sim.check_status()
if not status['output_exists']:
    print("Error: No output file found!")
    print(f"Expected: {sim.output_file}")
    sys.exit(1)

print(f"Analyzing simulation: {sim.run_dir}")
print(f"Current time: {status['current_time']:.2f} / {status['max_time']:.2f}")
print(f"Progress: {status['progress']:.1f}%")

# Load grids
print("\nLoading grids...")
x, kx, y, ky, z, p, j = sim.load_grids()
print(f"  x: {len(x)} points, range [{x[0]:.2f}, {x[-1]:.2f}]")
print(f"  y: {len(y)} points, range [{y[0]:.2f}, {y[-1]:.2f}]")
print(f"  kx: {len(kx)} modes")
print(f"  ky: {len(ky)} modes")

# Load time traces
print("\nLoading time traces...")
t, hflux = sim.load_time_trace('hflux_x')
t, pflux = sim.load_time_trace('pflux_x')
print(f"  Time points: {len(t)}")
print(f"  Time range: [{t[0]:.2f}, {t[-1]:.2f}]")

# Load 3D fields at specific time
print("\nLoading 3D fields at t=5.0...")
try:
    t3d, phi, tf = sim.load_field_3d('phi', time=5.0)
    print(f"  Loaded phi at t={tf:.2f}")
    print(f"  Shape: {phi.shape}")
    print(f"  Max |phi|: {np.abs(phi).max():.3e}")
    
    # Transform to real space
    phi_xy = sim.field_to_realspace(phi, iz=-1)
    print(f"  Real space shape: {phi_xy.shape}")
except Exception as e:
    print(f"  Error loading field: {e}")

# Create visualizations
print("\nCreating plots...")

# 1. Simple flux plot
fig1, axes1 = sim.plot_fluxes()
fig1.savefig('fluxes.png', dpi=150, bbox_inches='tight')
print("  Saved: fluxes.png")

# 2. Field slice
fig2, ax2, im2 = sim.plot_field_slice('phi', time=5.0)
fig2.savefig('phi_slice.png', dpi=150, bbox_inches='tight')
print("  Saved: phi_slice.png")

# 3. Summary plot
fig3, axes3 = sim.plot_summary(time=5.0)
fig3.savefig('summary.png', dpi=150, bbox_inches='tight')
print("  Saved: summary.png")

print("\nAnalysis complete!")
plt.show()
