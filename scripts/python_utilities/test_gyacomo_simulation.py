"""
Test script for GyacomoSimulation class

This script performs basic validation of the wrapper functionality
without actually running simulations.
"""

import sys
sys.path.append('..')
from gyacomo_simulation import GyacomoSimulation
import tempfile
import shutil
from pathlib import Path

def test_initialization():
    """Test basic initialization"""
    print("Test 1: Initialization...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        sim = GyacomoSimulation(run_dir=tmpdir, jobnum=0)
        
        assert sim.run_dir == Path(tmpdir)
        assert sim.jobnum == 0
        assert isinstance(sim.params, dict)
        assert 'BASIC' in sim.params
        assert 'GRID' in sim.params
        assert 'GEOMETRY' in sim.params
        
        print("  ✓ Initialization successful")

def test_parameter_modification():
    """Test parameter modification"""
    print("\nTest 2: Parameter modification...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        sim = GyacomoSimulation(run_dir=tmpdir)
        
        # Modify various parameters
        sim.params['BASIC']['tmax'] = 15.0
        sim.params['GRID']['Nx'] = 8
        sim.params['GEOMETRY']['delta'] = 0.5
        sim.params['SPECIES'][0]['k_T_'] = 7.0
        
        assert sim.params['BASIC']['tmax'] == 15.0
        assert sim.params['GRID']['Nx'] == 8
        assert sim.params['GEOMETRY']['delta'] == 0.5
        assert sim.params['SPECIES'][0]['k_T_'] == 7.0
        
        print("  ✓ Parameter modification successful")

def test_params_io():
    """Test reading and writing params files"""
    print("\nTest 3: Parameter file I/O...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        sim = GyacomoSimulation(run_dir=tmpdir)
        
        # Modify some parameters
        sim.params['BASIC']['tmax'] = 20.0
        sim.params['GEOMETRY']['q0'] = 3.0
        
        # Write to file
        sim.write_params_file()
        
        assert sim.params_file.exists()
        
        # Create new instance and load
        sim2 = GyacomoSimulation(run_dir=tmpdir, jobnum=0)
        
        # Verify loaded correctly
        assert sim2.params['BASIC']['tmax'] == 20.0
        assert sim2.params['GEOMETRY']['q0'] == 3.0
        
        print("  ✓ Parameter file I/O successful")

def test_directory_setup():
    """Test directory creation"""
    print("\nTest 4: Directory setup...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        run_dir = Path(tmpdir) / "test_sim"
        sim = GyacomoSimulation(run_dir=run_dir)
        
        assert not run_dir.exists()
        
        sim.setup_directory()
        
        assert run_dir.exists()
        assert (run_dir / "params.in").exists()
        
        print("  ✓ Directory setup successful")

def test_load_existing():
    """Test loading from existing simulation"""
    print("\nTest 5: Load existing simulation...")
    
    # Try to load the PT_rho0.9 simulation from attachments
    sim_dir = Path(__file__).parent.parent.parent.parent / \
              "simulations/check_instability_tcv_gkyl/PT_rho0.9"
    
    if not sim_dir.exists():
        print("  ⊘ Skipped (example simulation not found)")
        return
    
    sim = GyacomoSimulation(run_dir=sim_dir, jobnum=0)
    
    # Check that params were loaded
    assert sim.params['GEOMETRY']['q0'] == 2.8
    assert sim.params['GEOMETRY']['delta'] == 0.35
    
    # Check status
    status = sim.check_status()
    assert isinstance(status, dict)
    assert 'output_exists' in status
    
    print(f"  ✓ Loaded existing simulation")
    print(f"    Output exists: {status['output_exists']}")
    if status['current_time']:
        print(f"    Current time: {status['current_time']:.2f}")

def test_repr():
    """Test string representation"""
    print("\nTest 6: String representation...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        sim = GyacomoSimulation(run_dir=tmpdir, jobnum=0)
        repr_str = repr(sim)
        
        assert 'GyacomoSimulation' in repr_str
        assert tmpdir in repr_str
        assert 'job=0' in repr_str
        
        print(f"  ✓ {repr_str}")

def main():
    """Run all tests"""
    print("="*60)
    print("GyacomoSimulation Validation Tests")
    print("="*60)
    
    try:
        test_initialization()
        test_parameter_modification()
        test_params_io()
        test_directory_setup()
        test_load_existing()
        test_repr()
        
        print("\n" + "="*60)
        print("All tests passed! ✓")
        print("="*60)
        
    except AssertionError as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return 1
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
