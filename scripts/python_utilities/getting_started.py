#!/usr/bin/env python
"""
Getting Started with GyacomoSimulation

Interactive tutorial demonstrating the Python wrapper for Gyacomo.
Run this script to see basic usage patterns.
"""

from gyacomo_simulation import GyacomoSimulation
import numpy as np
import matplotlib.pyplot as plt

def main():
    print("="*70)
    print("GyacomoSimulation - Getting Started Tutorial")
    print("="*70)
    
    # =================================================================
    # PART 1: Create a new simulation
    # =================================================================
    print("\n" + "="*70)
    print("PART 1: Creating a New Simulation")
    print("="*70)
    
    print("\nCreating a simulation instance...")
    sim = GyacomoSimulation(
        run_dir='./tutorial_simulation',
        jobnum=0
    )
    print(f"✓ Created: {sim}")
    
    # =================================================================
    # PART 2: View default parameters
    # =================================================================
    print("\n" + "="*70)
    print("PART 2: Exploring Default Parameters")
    print("="*70)
    
    print("\nAvailable namelists:")
    for namelist in sim.params.keys():
        if namelist == 'SPECIES':
            print(f"  - {namelist} (list of {len(sim.params[namelist])} species)")
        else:
            n_params = len(sim.params[namelist])
            print(f"  - {namelist} ({n_params} parameters)")
    
    print("\nSample parameters:")
    print(f"  Time step:        {sim.params['BASIC']['dt']}")
    print(f"  Max time:         {sim.params['BASIC']['tmax']}")
    print(f"  Safety factor:    {sim.params['GEOMETRY']['q0']}")
    print(f"  Triangularity:    {sim.params['GEOMETRY']['delta']}")
    print(f"  Grid (Nx,Ny,Nz):  {sim.params['GRID']['Nx']}, {sim.params['GRID']['Ny']}, {sim.params['GRID']['Nz']}")
    
    print(f"\nIon parameters:")
    print(f"  k_N (density):    {sim.params['SPECIES'][0]['k_N_']}")
    print(f"  k_T (temperature): {sim.params['SPECIES'][0]['k_T_']}")
    
    # =================================================================
    # PART 3: Modify parameters
    # =================================================================
    print("\n" + "="*70)
    print("PART 3: Modifying Parameters")
    print("="*70)
    
    print("\nModifying simulation parameters...")
    
    # Basic parameters
    sim.params['BASIC']['tmax'] = 15.0
    sim.params['BASIC']['dt'] = 0.0005
    print(f"  ✓ Set tmax = {sim.params['BASIC']['tmax']}")
    
    # Grid resolution
    sim.params['GRID']['Nx'] = 4
    sim.params['GRID']['Ny'] = 32
    sim.params['GRID']['Nz'] = 32
    print(f"  ✓ Set grid = ({sim.params['GRID']['Nx']}, {sim.params['GRID']['Ny']}, {sim.params['GRID']['Nz']})")
    
    # Geometry (create ITG-unstable case)
    sim.params['GEOMETRY']['q0'] = 1.5
    sim.params['GEOMETRY']['delta'] = 0.0  # Circular
    print(f"  ✓ Set q0 = {sim.params['GEOMETRY']['q0']}, delta = {sim.params['GEOMETRY']['delta']}")
    
    # Species gradients (strong ITG drive)
    sim.params['SPECIES'][0]['k_T_'] = 9.0  # Strong ion gradient
    sim.params['SPECIES'][1]['k_T_'] = 2.0  # Weak electron gradient
    print(f"  ✓ Set ion k_T = {sim.params['SPECIES'][0]['k_T_']}")
    
    # Model settings
    sim.params['MODEL']['LINEARITY'] = 'linear'
    print(f"  ✓ Set mode = {sim.params['MODEL']['LINEARITY']}")
    
    # =================================================================
    # PART 4: Setup simulation directory
    # =================================================================
    print("\n" + "="*70)
    print("PART 4: Setting Up Simulation Directory")
    print("="*70)
    
    print(f"\nCreating directory: {sim.run_dir}")
    sim.setup_directory(clean=True)
    print(f"  ✓ Directory created")
    print(f"  ✓ Parameters written to: {sim.params_file}")
    
    # Show part of the params file
    print("\nFirst 20 lines of params.in:")
    print("-" * 70)
    with open(sim.params_file, 'r') as f:
        for i, line in enumerate(f):
            if i >= 20:
                print("  ...")
                break
            print(f"  {line.rstrip()}")
    print("-" * 70)
    
    # =================================================================
    # PART 5: How to run (without actually running)
    # =================================================================
    print("\n" + "="*70)
    print("PART 5: Execution (Demo - Not Actually Running)")
    print("="*70)
    
    print("\nTo run this simulation, you would call:")
    print(f"\n  sim.run(nproc=8, blocking=True)\n")
    print("This executes:")
    
    Nx_procs = sim.params['GRID']['Nx']
    Ny_procs = 8 // Nx_procs
    print(f"  mpirun -np 8 {sim.executable} {Nx_procs} {Ny_procs} 1 < params.in")
    
    print("\nFor this tutorial, we'll skip the actual execution.")
    print("See example_setup_run.py for a complete setup-and-run workflow.")
    
    # =================================================================
    # PART 6: Analysis (if output exists)
    # =================================================================
    print("\n" + "="*70)
    print("PART 6: Analysis Capabilities")
    print("="*70)
    
    print("\nOnce a simulation has run, you can analyze it with:")
    
    print("\n  # Check status")
    print("  status = sim.check_status()")
    print("  # Returns: {'output_exists', 'current_time', 'max_time', 'progress'}")
    
    print("\n  # Load grids")
    print("  x, kx, y, ky, z, p, j = sim.load_grids()")
    
    print("\n  # Load time traces")
    print("  t, hflux = sim.load_time_trace('hflux_x')")
    print("  t, pflux = sim.load_time_trace('pflux_x')")
    
    print("\n  # Load 3D fields")
    print("  t3d, phi, tf = sim.load_field_3d('phi', time=5.0)")
    
    print("\n  # Compute growth rates")
    print("  result = sim.compute_growth_rate(field_name='phi')")
    print("  gamma = np.real(result['omega'])")
    print("  freq = np.imag(result['omega'])")
    
    print("\n  # Create plots")
    print("  fig, axes = sim.plot_summary(time=5.0)")
    print("  fig, ax = sim.plot_dispersion()")
    
    # =================================================================
    # PART 7: Try analyzing existing simulation
    # =================================================================
    print("\n" + "="*70)
    print("PART 7: Analyzing Existing Simulation (if available)")
    print("="*70)
    
    # Try to find an existing simulation
    import os
    existing_sim = '/Users/ahoffman/gyacomo/simulations/check_instability_tcv_gkyl/PT_rho0.9'
    
    if os.path.exists(existing_sim):
        print(f"\nFound existing simulation: {existing_sim}")
        sim_existing = GyacomoSimulation(run_dir=existing_sim, jobnum=0)
        
        status = sim_existing.check_status()
        print(f"\nStatus:")
        print(f"  Output exists: {status['output_exists']}")
        
        if status['output_exists']:
            print(f"  Current time:  {status['current_time']:.2f}")
            print(f"  Max time:      {status['max_time']:.2f}")
            print(f"  Progress:      {status['progress']:.1f}%")
            
            # Load some data
            print("\nLoading grids...")
            x, kx, y, ky, z, p, j = sim_existing.load_grids()
            print(f"  Grid sizes: Nx={len(kx)}, Ny={len(ky)}, Nz={len(z)}")
            
            print("\nLoading time trace...")
            t, hflux = sim_existing.load_time_trace('hflux_x')
            print(f"  Time points: {len(t)}")
            print(f"  Time range: [{t[0]:.2f}, {t[-1]:.2f}]")
            
            print("\n✓ Successfully loaded existing simulation data!")
    else:
        print(f"\nNo existing simulation found at {existing_sim}")
        print("Run example_setup_run.py first to create a simulation.")
    
    # =================================================================
    # Summary
    # =================================================================
    print("\n" + "="*70)
    print("TUTORIAL COMPLETE")
    print("="*70)
    
    print("\nKey Takeaways:")
    print("  1. Create simulation: GyacomoSimulation(run_dir, jobnum)")
    print("  2. Modify parameters: sim.params[namelist][parameter] = value")
    print("  3. Setup directory: sim.setup_directory()")
    print("  4. Run simulation: sim.run(nproc=N)")
    print("  5. Analyze results: sim.load_*() and sim.plot_*() methods")
    
    print("\nNext Steps:")
    print("  • See examples/ directory for complete workflows")
    print("  • Read README_GyacomoSimulation.md for full documentation")
    print("  • Try example_setup_run.py to run a real simulation")
    print("  • Use example_analysis.py to analyze existing results")
    
    print("\n" + "="*70)

if __name__ == '__main__':
    main()
