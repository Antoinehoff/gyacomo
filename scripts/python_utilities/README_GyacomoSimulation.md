# GyacomoSimulation - Python Wrapper for Gyacomo

A comprehensive Python interface for managing Gyacomo gyrokinetic turbulence simulations. This module provides tools to configure parameters, execute simulations via MPI, and analyze HDF5 output data.

## Features

- 📝 **Parameter Management**: Read/write Fortran namelist files with Python dictionaries
- 🚀 **MPI Execution**: Run simulations with automatic process decomposition
- 📊 **Data Analysis**: Load time traces, 3D fields, and grids from HDF5 outputs
- 📈 **Visualization**: Built-in plotting for fluxes, fields, and dispersion relations
- 🔄 **Restart Capability**: Continue simulations from previous runs
- 📦 **Batch Processing**: Set up parameter scans and multi-case studies

## Quick Start

```python
from gyacomo_simulation import GyacomoSimulation

# Create a new simulation
sim = GyacomoSimulation(run_dir='./my_simulation', jobnum=0)

# Modify parameters
sim.params['BASIC']['tmax'] = 10.0
sim.params['GEOMETRY']['delta'] = 0.3  # Set triangularity
sim.params['SPECIES'][0]['k_T_'] = 6.0 # Ion temperature gradient

# Setup and run
sim.setup_directory()
sim.run(nproc=8)

# Analyze results
fig, axes = sim.plot_summary(time=5.0)
result = sim.compute_growth_rate(field_name='phi')
```

## Installation

No installation required - just ensure the Python utilities are in your path:

```bash
cd /path/to/gyacomo/scripts/python_utilities
python
>>> from gyacomo_simulation import GyacomoSimulation
```

### Dependencies

Required Python packages:
- `numpy` - Array operations
- `matplotlib` - Visualization
- `h5py` - HDF5 file I/O (used by load_data.py)

Existing Gyacomo utilities (in same directory):
- `load_data.py` - HDF5 data loading
- `tools.py` - Analysis utilities
- `fourier.py` - FFT operations

## Class Overview

### `GyacomoSimulation`

Main class for managing simulations.

#### Initialization

```python
sim = GyacomoSimulation(run_dir, jobnum=0, executable=None)
```

**Parameters:**
- `run_dir`: Directory for simulation files
- `jobnum`: Job index for restarts (default: 0)
- `executable`: Path to gyacomo binary (auto-detected if None)

**Attributes:**
- `params`: Dictionary of all namelist parameters
- `run_dir`: Simulation directory path
- `output_file`: Path to outputs_XX.h5
- `executable`: Path to gyacomo executable

#### Parameter Structure

The `params` dictionary mirrors Fortran namelists:

```python
sim.params = {
    'BASIC': {'nrun', 'dt', 'tmax', 'maxruntime', 'job2load'},
    'GRID': {'pmax', 'jmax', 'Nx', 'Lx', 'Ny', 'Ly', 'Nz', ...},
    'GEOMETRY': {'geom', 'q0', 'shear', 'eps', 'kappa', 'delta', ...},
    'DIAGNOSTICS': {'dtsave_0d', 'dtsave_3d', 'dtsave_5d', ...},
    'MODEL': {'LINEARITY', 'Na', 'nu', 'beta', 'ADIAB_E', ...},
    'CLOSURE': {'hierarchy_closure', 'dmax', 'nonlinear_closure', ...},
    'SPECIES': [{name_, tau_, sigma_, q_, k_N_, k_T_}, ...],
    'COLLISION': {'collision_model', 'GK_CO'},
    'INITIAL': {'INIT_OPT'},
    'TIME_INTEGRATION': {'numerical_scheme'},
    'UNITS': {'n_ref', 'T_ref', 'R_ref', 'B_ref', ...}
}
```

## Methods Reference

### Setup and Execution

#### `setup_directory(clean=False)`
Create run directory and write params.in file.

```python
sim.setup_directory(clean=True)  # Remove existing directory
```

#### `run(nproc=8, blocking=True, verbose=True)`
Execute simulation with MPI.

```python
result = sim.run(nproc=8, blocking=True)  # Wait for completion
process = sim.run(nproc=8, blocking=False)  # Run in background
```

**Parameters:**
- `nproc`: Number of MPI processes
- `blocking`: Wait for completion
- `verbose`: Print progress

#### `check_status()`
Check simulation progress.

```python
status = sim.check_status()
# Returns: {'output_exists', 'current_time', 'max_time', 'progress'}
```

#### `write_params_file(filename=None)`
Write current parameters to namelist file.

```python
sim.write_params_file('custom_params.in')
```

#### `load_params_from_file(filename)`
Load parameters from existing namelist file.

```python
sim.load_params_from_file('params.in')
```

### Data Loading

#### `load_params()`
Load parameters from HDF5 output.

```python
params = sim.load_params()
```

#### `load_grids()`
Load coordinate grids.

```python
x, kx, y, ky, z, p, j = sim.load_grids()
```

#### `load_time_trace(varname)`
Load 0D time series data.

```python
t, hflux = sim.load_time_trace('hflux_x')
t, pflux = sim.load_time_trace('pflux_x')
```

Available variables: `'hflux_x'`, `'pflux_x'`, `'gflux_x'`, `'energy'`

#### `load_field_3d(varname, time=None)`
Load 3D spectral field at specific time.

```python
t3d, phi, tf = sim.load_field_3d('phi', time=5.0)
```

Available fields: `'phi'`, `'psi'`, `'Na00'`, `'dens'`, `'upar'`, `'Tpar'`, `'temp'`

#### `field_to_realspace(field_kxkyz, iz=-1)`
Transform spectral field to real space.

```python
phi_xy = sim.field_to_realspace(phi_3d, iz=-1)  # Outboard midplane
```

### Analysis

#### `compute_growth_rate(field_name='phi', kx_val=0.0, ky_index=None, time_range=None, iz=None)`
Compute complex growth rate from field evolution.

```python
result = sim.compute_growth_rate(
    field_name='phi',
    kx_val=0.0,
    time_range=(2.0, 5.0)
)
# Returns: {'ky': array, 'omega': array (complex), 'time': array}

gamma = np.real(result['omega'])  # Growth rate
omega = np.imag(result['omega'])  # Frequency
```

### Visualization

#### `plot_fluxes(ax=None, species_labels=None)`
Plot heat and particle flux time traces.

```python
fig, axes = sim.plot_fluxes()
```

#### `plot_field_slice(field_name, time=None, iz=-1, ax=None, cmap='seismic', **kwargs)`
Plot 2D real-space field slice.

```python
fig, ax, im = sim.plot_field_slice('phi', time=5.0, iz=-1)
```

#### `plot_summary(time=None, figsize=(12, 8))`
Create 4-panel summary plot (fluxes + fields).

```python
fig, axes = sim.plot_summary(time=5.0)
```

#### `plot_dispersion(field_name='phi', kx_val=0.0, ky_range=None, time_range=None, figsize=(10, 6))`
Plot growth rate and frequency vs wavenumber.

```python
fig, ax = sim.plot_dispersion(
    field_name='phi',
    ky_range=(0.1, 1.0),
    time_range=(2.0, 5.0)
)
```

## Examples

See the [examples/](examples/) directory for complete scripts:

1. **[example_setup_run.py](examples/example_setup_run.py)** - Setup and execute simulations
2. **[example_analysis.py](examples/example_analysis.py)** - Load and analyze results
3. **[example_growth_rate_scan.py](examples/example_growth_rate_scan.py)** - Compute dispersion relations
4. **[example_parameter_scan.py](examples/example_parameter_scan.py)** - Parameter scan setup

## Common Workflows

### 1. Linear ITG Simulation

```python
# Setup ITG-unstable case
sim = GyacomoSimulation(run_dir='./itg_test')
sim.params['MODEL']['LINEARITY'] = 'linear'
sim.params['SPECIES'][0]['k_T_'] = 6.0  # Strong ion gradient
sim.params['SPECIES'][1]['k_T_'] = 2.0  # Weak electron gradient
sim.setup_directory()
sim.run(nproc=8)

# Analyze growth rate
result = sim.compute_growth_rate()
print(f"Max gamma: {np.real(result['omega']).max():.4f}")
```

### 2. Triangularity Comparison

```python
# Compare positive and negative triangularity
for delta in [0.3, -0.3]:
    run_dir = f'./tri_scan/delta_{delta:.1f}'
    sim = GyacomoSimulation(run_dir=run_dir)
    sim.params['GEOMETRY']['delta'] = delta
    sim.setup_directory()
    sim.run(nproc=8)
```

### 3. Restart Long Simulation

```python
# Initial run (0-10 time units)
sim = GyacomoSimulation(run_dir='./long_run', jobnum=0)
sim.params['BASIC']['tmax'] = 10.0
sim.run(nproc=8)

# Restart (10-20 time units)
sim_restart = GyacomoSimulation(run_dir='./long_run', jobnum=1)
sim_restart.params['BASIC']['job2load'] = 0
sim_restart.params['BASIC']['tmax'] = 20.0
sim_restart.run(nproc=8)
```

### 4. Batch Analysis

```python
import glob
import numpy as np

# Analyze all simulations in directory
results = []
for sim_dir in glob.glob('./parameter_scan/*/'):
    sim = GyacomoSimulation(run_dir=sim_dir)
    status = sim.check_status()
    
    if status['output_exists']:
        result = sim.compute_growth_rate()
        gamma_max = np.real(result['omega']).max()
        results.append({'dir': sim_dir, 'gamma': gamma_max})

# Print summary
for r in sorted(results, key=lambda x: x['gamma'], reverse=True):
    print(f"{r['dir']}: gamma = {r['gamma']:.4f}")
```

## Parameter Guide

### Critical Physics Parameters

**Density gradient** (`k_N_`): Affects drive for instabilities
```python
sim.params['SPECIES'][0]['k_N_'] = 3.0  # R/L_n
```

**Temperature gradient** (`k_T_`): Primary ITG drive
```python
sim.params['SPECIES'][0]['k_T_'] = 6.0  # R/L_T
```

**Safety factor** (`q0`): Controls mode structure
```python
sim.params['GEOMETRY']['q0'] = 2.0
```

**Magnetic shear** (`shear`): Affects radial mode structure
```python
sim.params['GEOMETRY']['shear'] = 1.0  # s = r/q * dq/dr
```

**Triangularity** (`delta`): Geometry shape parameter
```python
sim.params['GEOMETRY']['delta'] = 0.3  # Positive triangularity
```

**Collision frequency** (`nu`): Damping parameter
```python
sim.params['MODEL']['nu'] = 0.1  # Normalized to c_s/R
```

### Resolution Parameters

**Radial modes** (`Nx`): Perpendicular turbulence structure
```python
sim.params['GRID']['Nx'] = 4  # Linear: 2-4, Nonlinear: 32-128
```

**Binormal modes** (`Ny`): Perpendicular turbulence structure
```python
sim.params['GRID']['Ny'] = 32  # Linear: 24-32, Nonlinear: 128-256
```

**Parallel points** (`Nz`): Parallel resolution
```python
sim.params['GRID']['Nz'] = 32  # Typically 24-64
```

**Velocity space** (`pmax`, `jmax`): Distribution function resolution
```python
sim.params['GRID']['pmax'] = 4  # Hermite degree (parallel)
sim.params['GRID']['jmax'] = 2  # Laguerre degree (perpendicular)
```

## Troubleshooting

### Common Issues

**"No output file found"**
- Simulation hasn't run yet or failed
- Check `sim.check_status()` and stdout file

**MPI decomposition errors**
- Ensure `nproc = Nx_procs * Ny_procs * 1`
- Adjust Nx parameter: `nproc // Ny` must be integer

**Import errors**
- Ensure load_data.py, tools.py, fourier.py are in same directory
- Check Python path includes utilities directory

**HDF5 errors**
- Verify h5py is installed: `pip install h5py`
- Check file permissions on output files

**Growth rate calculation fails**
- Not enough time frames saved (increase diagnostic frequency)
- Field doesn't exist (check available fields)

### Debug Tips

```python
# Check what's available
print(sim.params)
print(f"Output file: {sim.output_file}")
print(f"Exists: {sim.output_file.exists()}")

# Inspect HDF5 structure
import h5py
with h5py.File(sim.output_file, 'r') as f:
    print(list(f['data'].keys()))
```

## Advanced Usage

### Custom Parameter Initialization

```python
# Start from existing simulation
sim_new = GyacomoSimulation(run_dir='./new_run')
sim_new.load_params_from_file('../old_run/params.in')
sim_new.params['GEOMETRY']['delta'] = 0.5  # Modify one parameter
sim_new.setup_directory()
```

### Non-blocking Execution

```python
import time

# Start simulation in background
process = sim.run(nproc=8, blocking=False)

# Monitor progress
while process.poll() is None:
    time.sleep(10)
    status = sim.check_status()
    if status['progress']:
        print(f"Progress: {status['progress']:.1f}%")
```

### Custom Analysis Pipeline

```python
class MyAnalysis(GyacomoSimulation):
    def compute_custom_metric(self):
        """Add custom analysis method"""
        t, hflux = self.load_time_trace('hflux_x')
        return np.trapz(hflux, t)  # Integrated flux
    
    def plot_custom(self):
        """Add custom plot"""
        # Your plotting code here
        pass
```

## Contributing

To extend functionality:
1. Add methods to `GyacomoSimulation` class
2. Create example scripts in `examples/`
3. Update this README
4. Test with various simulation types

## License

This wrapper follows the same license as Gyacomo. See main repository for details.

## Contact

For issues specific to this Python wrapper, check:
- Main Gyacomo documentation: `../../README.md`
- Example scripts: `examples/`
- Existing analysis notebooks: `analysis.ipynb`

---

**Version:** 1.0  
**Last updated:** February 2026  
**Compatible with:** Gyacomo 2023
