# Gyacomo Python Wrapper - Examples

This directory contains example scripts demonstrating how to use the `GyacomoSimulation` class for setting up, running, and analyzing Gyacomo gyrokinetic simulations.

## Prerequisites

Make sure the following modules are in your Python path:
- `gyacomo_simulation.py` (the main wrapper class)
- `load_data.py` (HDF5 data loading utilities)
- `tools.py` (Fourier transforms and analysis tools)
- `fourier.py` (FFT utilities)

## Examples

### 1. `example_setup_run.py`
**Purpose:** Create and configure a new simulation from scratch

**What it demonstrates:**
- Creating a `GyacomoSimulation` instance
- Modifying parameters (grid, geometry, species, etc.)
- Setting up the run directory
- Executing the simulation with MPI
- Checking simulation status

**Usage:**
```bash
cd examples
python example_setup_run.py
```

**Key features:**
- Shows how to set ITG-like parameters
- Demonstrates parameter modification
- Includes safety comments (run command is commented out)

---

### 2. `example_analysis.py`
**Purpose:** Load and analyze existing simulation results

**What it demonstrates:**
- Loading an existing simulation
- Extracting time traces (heat flux, particle flux)
- Loading 3D field data
- Creating basic visualizations

**Usage:**
```bash
python example_analysis.py
```

**Outputs:**
- `fluxes.png` - Heat and particle flux time traces
- `phi_slice.png` - Electrostatic potential at outboard midplane
- `summary.png` - 4-panel overview plot

**Customize:**
Change the `run_dir` parameter to point to your simulation directory.

---

### 3. `example_growth_rate_scan.py`
**Purpose:** Compute and compare growth rates (replicates analysis.ipynb workflow)

**What it demonstrates:**
- Computing complex growth rates from field evolution
- Comparing multiple simulations (PT vs NT triangularity)
- Creating publication-quality dispersion plots
- Extracting peak growth rates and frequencies

**Usage:**
```bash
python example_growth_rate_scan.py
```

**Outputs:**
- `growth_rate_comparison_rho0.9.png` - Dispersion relation comparison
- Console summary of max growth rates

**Features:**
- Analyzes phi field at kx=0 for all ky modes
- Automatically determines appropriate time window
- Formats plots for publication

---

### 4. `example_parameter_scan.py`
**Purpose:** Set up multiple simulations for parameter scans

**What it demonstrates:**
- Creating multiple simulation directories
- Scanning over a parameter (density gradient k_N)
- Organizing results systematically
- Batch analysis template

**Usage:**
```bash
python example_parameter_scan.py
```

**What it does:**
- Creates 5 simulation directories with different k_N values
- Writes params.in file for each case
- Provides instructions for running and analyzing

**Customize:**
- Change `k_N_values` to scan different parameters
- Modify base parameters to match your physics case
- Uncomment run loop for sequential execution

---

## Quick Start

### Basic Workflow

1. **Setup a simulation:**
```python
from gyacomo_simulation import GyacomoSimulation

sim = GyacomoSimulation(run_dir='./my_simulation', jobnum=0)
sim.params['BASIC']['tmax'] = 10.0
sim.params['GRID']['Nx'] = 4
sim.setup_directory()
```

2. **Run the simulation:**
```python
sim.run(nproc=8, blocking=True)
```

3. **Analyze results:**
```python
# Check status
status = sim.check_status()
print(f"Progress: {status['progress']:.1f}%")

# Load data
t, hflux = sim.load_time_trace('hflux_x')
x, kx, y, ky, z, p, j = sim.load_grids()

# Visualize
fig, axes = sim.plot_summary(time=5.0)
```

4. **Compute growth rates:**
```python
result = sim.compute_growth_rate(field_name='phi', kx_val=0.0)
gamma = np.real(result['omega'])
print(f"Max growth rate: {gamma.max():.4f}")
```

## Common Tasks

### Modify geometry parameters
```python
sim.params['GEOMETRY']['q0'] = 2.0      # Safety factor
sim.params['GEOMETRY']['shear'] = 1.5   # Magnetic shear
sim.params['GEOMETRY']['delta'] = 0.3   # Triangularity
sim.params['GEOMETRY']['kappa'] = 1.6   # Elongation
```

### Change species gradients
```python
sim.params['SPECIES'][0]['k_N_'] = 3.0  # Ion density gradient
sim.params['SPECIES'][0]['k_T_'] = 6.0  # Ion temperature gradient
sim.params['SPECIES'][1]['k_N_'] = 3.0  # Electron density gradient
sim.params['SPECIES'][1]['k_T_'] = 3.0  # Electron temperature gradient
```

### Adjust resolution
```python
sim.params['GRID']['Nx'] = 4     # Radial Fourier modes
sim.params['GRID']['Ny'] = 32    # Binormal Fourier modes
sim.params['GRID']['Nz'] = 32    # Parallel grid points
sim.params['GRID']['pmax'] = 6   # Hermite polynomial degree
sim.params['GRID']['jmax'] = 3   # Laguerre polynomial degree
```

### Configure diagnostics
```python
sim.params['DIAGNOSTICS']['dtsave_0d'] = 0.05  # Time traces
sim.params['DIAGNOSTICS']['dtsave_3d'] = 0.2   # 3D fields
sim.params['DIAGNOSTICS']['dtsave_5d'] = 1.0   # Full distribution
```

## Troubleshooting

**Q: Import errors for load_data or tools**
- Ensure you're running from the examples directory
- Check that `sys.path.append('..')` points to python_utilities

**Q: Executable not found**
- Set explicitly: `sim.executable = '/path/to/gyacomo23_dp'`
- Or ensure executables are in `../../bin/`

**Q: No output file found**
- Check that simulation has run: `sim.check_status()`
- Verify output file location: `print(sim.output_file)`

**Q: MPI decomposition errors**
- Ensure `nproc = Nx_procs * Ny_procs`
- Adjust Nx and Ny parameters accordingly

## Advanced Usage

### Restart simulations
```python
# Initial run
sim = GyacomoSimulation(run_dir='./my_sim', jobnum=0)
sim.run(nproc=8)

# Restart from previous job
sim_restart = GyacomoSimulation(run_dir='./my_sim', jobnum=1)
sim_restart.params['BASIC']['job2load'] = 0  # Load from job 0
sim_restart.params['BASIC']['tmax'] = 20.0   # Run longer
sim_restart.run(nproc=8)
```

### Custom plotting
```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
t, flux = sim.load_time_trace('hflux_x')
ax.semilogy(t, np.abs(flux))
ax.set_xlabel('Time')
ax.set_ylabel('Heat Flux')
plt.savefig('custom_plot.png')
```

### Batch processing
```python
# Analyze multiple simulations
results = []
for case_name in ['case1', 'case2', 'case3']:
    sim = GyacomoSimulation(run_dir=f'./sims/{case_name}')
    result = sim.compute_growth_rate()
    results.append({'case': case_name, 'gamma': np.real(result['omega']).max()})
```

## Additional Resources

- Main class documentation: `../gyacomo_simulation.py`
- Data loading utilities: `../load_data.py`
- Analysis tools: `../tools.py`
- Gyacomo documentation: `../../../README.md`

## Contributing

To add new examples:
1. Create a new file `example_<name>.py`
2. Include docstring explaining purpose
3. Add comments for clarity
4. Update this README

---

**Last updated:** February 2026
