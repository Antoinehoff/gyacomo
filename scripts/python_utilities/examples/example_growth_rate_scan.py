"""
Example: Growth Rate Scan with Triangularity Comparison

This script demonstrates how to:
1. Analyze growth rates from simulations
2. Compare different cases (PT vs NT triangularity)
3. Create publication-quality dispersion plots

This replicates the analysis from analysis.ipynb cell 2.
"""

import sys
sys.path.append('..')
from gyacomo_simulation import GyacomoSimulation
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "text.usetex": False,
    "font.size": 14
})

# Define cases to compare
rho = 0.9
tri_list = ['PT', 'NT']

# Storage for results
omega_dict = {}

# Analyze each triangularity case
for tri_type in tri_list:
    print(f"\nAnalyzing {tri_type} triangularity at rho={rho}...")
    
    # Load simulation
    run_dir = f'../../../simulations/check_instability_tcv_gkyl/{tri_type}_rho{rho}'
    sim = GyacomoSimulation(run_dir=run_dir, jobnum=0)
    
    # Check if output exists
    status = sim.check_status()
    if not status['output_exists']:
        print(f"  Warning: No output found for {tri_type}")
        continue
    
    # Load grids
    x, kx, y, ky, z, p, j = sim.load_grids()
    print(f"  Grid: Nz={len(z)}, Nkx={len(kx)}, Nky={len(ky)}")
    
    # Compute growth rates for kx=0, all ky
    # Use time range from t/4 to t/2 for fitting
    result = sim.compute_growth_rate(
        field_name='phi',
        kx_val=0.0,
        time_range=None  # Auto-determines from available data
    )
    
    omega_dict[tri_type] = result['omega']
    
    print(f"  Computed growth rates for {len(result['ky'])} modes")
    print(f"  Max gamma: {np.real(result['omega']).max():.3f}")
    print(f"  ky at max gamma: {result['ky'][np.argmax(np.real(result['omega']))]:.3f}")

# Plot comparison
print("\nCreating comparison plot...")
fig, ax = plt.subplots(figsize=(8, 6))

ik0 = 0  # Start from first mode (can skip ky=0 if needed)

for tri_type in tri_list:
    if tri_type not in omega_dict:
        continue
    
    # Reload to get ky values (they should be same for all cases)
    run_dir = f'../../../simulations/check_instability_tcv_gkyl/{tri_type}_rho{rho}'
    sim = GyacomoSimulation(run_dir=run_dir, jobnum=0)
    x, kx, y, ky, z, p, j = sim.load_grids()
    
    omega = omega_dict[tri_type]
    
    # Plot growth rate
    ax.plot(ky[ik0:], np.real(omega[ik0:]), 
           label=f'{tri_type} $\gamma$', linestyle='-', marker='o')
    
    # Plot frequency (scaled for visibility)
    ax.plot(ky[ik0:], np.imag(omega[ik0:])/4, 
           label=f'{tri_type} $\omega/4$', linestyle='--', marker='s')

ax.set_xlabel(r'$k_y \rho_s$')
ax.set_ylabel(r'$c_s/R$')
ax.set_title(f'Growth Rate Comparison at $\\rho={rho}$')
ax.set_xlim(0.1, 1.0)
ax.set_ylim(-5, 7)
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'growth_rate_comparison_rho{rho}.png', dpi=300, bbox_inches='tight')
print(f"Saved: growth_rate_comparison_rho{rho}.png")

plt.show()

# Print summary
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
for tri_type in tri_list:
    if tri_type not in omega_dict:
        continue
    omega = omega_dict[tri_type]
    gamma = np.real(omega)
    freq = np.imag(omega)
    
    print(f"\n{tri_type} Triangularity:")
    print(f"  Max growth rate: {gamma.max():.4f} c_s/R")
    print(f"  ky at max growth: {ky[np.argmax(gamma)]:.4f} rho_s")
    print(f"  Frequency range: [{freq.min():.4f}, {freq.max():.4f}] c_s/R")
