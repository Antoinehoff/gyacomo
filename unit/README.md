# Unit tests

Run them with

```sh
make check          # from the repository root
```

Each `test_*.f90` here is a standalone program that exercises one module and
exits non-zero if any assertion fails. `make check` builds them all and reports
a summary; the exit status is non-zero if any suite failed.

If your MPI needs a launcher even for a single rank:

```sh
make check TESTRUNNER="mpirun -np 1"
```

## What is covered

| suite | module | what it pins down |
|---|---|---|
| `test_prec_const` | `prec_const_mod` | the hand-typed constants against recomputed values |
| `test_utility` | `utility_mod` | `decomp1D` partitions exhaustively; `is_nan` / `is_inf` |
| `test_lag_interp` | `lag_interp_mod` | cubic exactness of the interpolator and its derivative, on ascending and descending grids, real, complex and 2D |
| `test_calculus` | `calculus_mod` | polynomial exactness and order of accuracy of every z stencil, the staggered half-cell offsets, and the too-few-ghosts guard |
| `test_numerics` | `numerics_mod` | `taylor_kernel_n` against the closed form `y^n exp(-y)/n!` |
| `test_coeff` | `coeff_mod` | Laguerre triple products `d_nks` against exact rational values, plus the symmetry that `build_dnjs_table` assumes |
| `test_grid` | `grid_mod` | spacings, Fourier mode ordering, ghost cell layout, index shortcuts, Simpson weights and `simpson_rule_z` |
| `test_time_integration` | `time_integration_mod` | every Butcher tableau against the Runge-Kutta order conditions for the order its name claims |
| `test_ExB_shear_flow` | `ExB_shear_flow_mod` | the shear setup, the kx remap of the fields against a reference built from the wavenumbers, and the shift that `Restore_ExB_shear_flow` rebuilds on a restart |

## How the build works

`make check` re-invokes make with `OBJDIR` and `MODDIR` pointing at
`obj/unit` and `mod/unit`, so the test build never shares state with whatever
`make dp` or `make debug` last left in `obj/` and `mod/`. Runtime checks
(`-fcheck=all`) are on, so an out of bounds access in the code under test fails
the test instead of going unnoticed.

Each test links against every object of the solver except `main.o`. That means
a test may `USE` any module without the Makefile needing a per-test dependency
list.

## Adding a test

1. Write `unit/test_<something>.f90` as a `PROGRAM` that uses `test_harness`:

   ```fortran
   PROGRAM test_something
     USE prec_const, ONLY: xp
     USE test_harness
     IMPLICIT NONE
     CALL start_suite('something')
     CALL section('what this group checks')
     CALL check_real('a value', got, expected)
     CALL finish_suite
   END PROGRAM test_something
   ```

2. Add `test_<something>` to `UNIT_PROGS` in the Makefile.

The harness offers `check`, `check_int`, `check_real`, `check_real_arr`,
`check_cmpx`, `check_cmpx_arr`, and `tol_of(factor)` which returns
`factor*EPSILON(1._xp)` so tolerances stay meaningful in the single precision
build too. Tolerances default to `tol_of(1.0e3_xp)` and are relative, measured
against `max(1, |expected|)`.

If your test needs the Gyacomo communicators (`my_id`, `num_procs_z`, ...),
call `ppinit` *before* `start_suite`; the harness only initializes MPI if
nobody else has.

To feed a namelist to a `*_readinputs` routine, open `lu_in` as a scratch file
and write the namelist into it rather than leaving a file on disk — see
`test_grid.f90`.

## A caveat worth knowing

`allocate_array` in `basic_mod` allocates without deallocating first, so a
`*_readinputs` routine that allocates can only be called once per process.
That is why `test_time_integration` takes the scheme on the command line and
`make check` runs it once per scheme instead of looping inside the program.
`test_ExB_shear_flow` takes its linearity the same way, for the same reason
plus `Setup_ExB_shear_flow` allocating its own arrays unconditionally; the
modes it is run with are in `EXB_MODES` in the Makefile.

That test also stands up a whole run, `readinputs` then `auxval`, because the
shear touches the grids and the operators. It does *not* call `basic_data`:
that goes looking for an input file on disk and would read the command line
argument as a file index. `dt` therefore keeps the default of `basic_mod`, and
the test picks its shearing rate against whatever that default is.
