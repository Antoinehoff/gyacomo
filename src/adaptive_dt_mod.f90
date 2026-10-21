MODULE adaptive_dt
  ! Adaptive time stepping at a chosen fraction of the stability limit of the
  ! explicit Runge-Kutta scheme in use.
  !
  ! The time step is bounded by two things:
  !  1. the linear operator of the run (parallel streaming, drifts, collisions,
  !     hyperdiffusion, FLR and field couplings), whose extreme eigenvalues are
  !     measured once at startup by a matrix-free Arnoldi iteration on the very
  !     RHS the time stepper uses (assemble_RHS without the Poisson bracket).
  !     The operator is block diagonal in (kx,ky) (or in ky when the parallel
  !     boundary condition chains kx modes), so every block gets its own Krylov
  !     space and all of them advance together within one RHS evaluation;
  !  2. the ExB advection, re-measured after every step from the real-space
  !     gradients of the advecting potentials, as a CFL frequency
  !     omega_NL = max_r ( kx_eff |d_y F| + ky_eff |d_x F| ).
  ! Both are combined through the stability region of the scheme itself: the
  ! Ritz spectrum shifted by +/- i omega_NL is tested against |R(h lambda)| <= 1,
  ! where R is the stability polynomial built from the Butcher tableau. The
  ! result is tabulated once in omega_NL so that each step only interpolates.
  ! The step is then dt = cfl_frac * h_stab(omega_NL), with gentle growth,
  ! immediate decrease, bounds and an exact landing on tmax.
  USE prec_const, ONLY: xp, imagu, mpi_xp_r, mpi_xp_c
  IMPLICIT NONE
  PRIVATE

  !! Namelist ADAPTIVE_DT (optional group, position independent)
  LOGICAL,  PUBLIC, PROTECTED :: use_adaptive_dt   = .FALSE.   ! master switch
  REAL(xp), PUBLIC, PROTECTED :: cfl_frac          = 1.0_xp    ! fraction of the stability limit to run at
  REAL(xp), PUBLIC, PROTECTED :: dt_min            = 1.0e-6_xp ! the run stops below this
  REAL(xp), PUBLIC, PROTECTED :: dt_max            = 1.0e+6_xp ! the step never exceeds this
  REAL(xp), PUBLIC, PROTECTED :: dt_growth         = 1.2_xp    ! largest growth factor per step
  REAL(xp), PUBLIC, PROTECTED :: dt_hyst           = 0.1_xp    ! no growth unless the target exceeds dt by this fraction
  LOGICAL,  PUBLIC, PROTECTED :: estimate_dt_lin   = .TRUE.    ! measure the linear limit (Arnoldi) at startup
  REAL(xp), PUBLIC, PROTECTED :: dt_lin            = -1.0_xp   ! linear limit; input when not estimated, output otherwise
  REAL(xp), PUBLIC, PROTECTED :: dt_lin_safety     = 0.9_xp    ! margin on the measured linear spectrum
  INTEGER,  PUBLIC, PROTECTED :: krylov_dim        = 24        ! Arnoldi dimension per block
  REAL(xp), PUBLIC, PROTECTED :: krylov_mem_levels = 4.0_xp    ! Krylov basis memory budget, in moments time levels
  INTEGER,  PUBLIC, PROTECTED :: n_omega_tab       = 32        ! size of the h_stab(omega) table
  LOGICAL,  PUBLIC, PROTECTED :: dt_verbose        = .FALSE.   ! print every dt change

  !! State, readable by the diagnostics
  REAL(xp), PUBLIC, PROTECTED :: omega_nl = 0._xp  ! ExB CFL frequency of the current state
  REAL(xp), PUBLIC, PROTECTED :: dt_est   = 0._xp  ! stability limit for the current state (cfl_frac = 1)
  REAL(xp), PUBLIC, PROTECTED :: dt_used  = 0._xp  ! dt of the last completed step
  REAL(xp), PUBLIC, PROTECTED :: S_real   = 0._xp  ! stability interval of the scheme on the real axis
  REAL(xp), PUBLIC, PROTECTED :: S_imag   = 0._xp  ! stability interval of the scheme on the imaginary axis
  INTEGER,  PUBLIC, PROTECTED :: n_dt_changes     = 0
  INTEGER,  PUBLIC, PROTECTED :: n_cfl_violations = 0      ! steps that exceeded the limit of their end state by > 1 %
  INTEGER,  PUBLIC, PROTECTED :: n_cfl_warnings   = 0      ! ... by > 10 %
  REAL(xp), PUBLIC, PROTECTED :: max_cfl_excess   = 0._xp  ! largest dt/limit - 1 seen
  INTEGER,  PUBLIC, PROTECTED :: n_rhs_evals      = 0  ! RHS evaluations spent by the Arnoldi
  REAL(xp), PUBLIC, PROTECTED :: arnoldi_time     = 0._xp
  ! The block that limits the linear time step and its eigenvalue
  REAL(xp),    PUBLIC, PROTECTED :: lim_kx = 0._xp, lim_ky = 0._xp, lim_resid = 0._xp
  COMPLEX(xp), PUBLIC, PROTECTED :: lim_eig = (0._xp, 0._xp)
  LOGICAL,     PUBLIC, PROTECTED :: row_blocks = .FALSE. ! blocks are ky rows rather than (kx,ky) modes
  LOGICAL,     PUBLIC, PROTECTED :: adt_stop   = .FALSE. ! the stable step fell below dt_min (checked by tesend)
  ! Linear limit per local (ky,kx) mode, for the output file
  REAL(xp), DIMENSION(:,:), ALLOCATABLE, PUBLIC, PROTECTED :: dt_lin_map

  !! Private state
  REAL(xp), PARAMETER :: HUGE_H   = 1.0e30_xp ! "no limit"
  REAL(xp), PARAMETER :: RITZ_TOL = 0.1_xp    ! relative residual a Ritz value must reach to count
  INTEGER,  PARAMETER :: N_KEEP   = 6         ! extreme Ritz values kept per block for the omega table
  REAL(xp) :: R_max = 0._xp                   ! radius of the smallest disk containing the stability region
  REAL(xp) :: flutter_factor = 0._xp          ! velocity bound of the magnetic flutter bracket
  REAL(xp) :: omega_tab_max = 0._xp, domega_tab = 0._xp
  REAL(xp), DIMENSION(:), ALLOCATABLE :: inv_h_tab    ! 1/h_stab on the omega grid (0 = no limit)
  COMPLEX(xp), DIMENSION(:), ALLOCATABLE :: ritz_list ! kept (inflated) Ritz values of the local blocks
  INTEGER :: n_ritz = 0
  LOGICAL :: dt_lin_estimated = .FALSE.
  ! Saved state around the Arnoldi pass
  COMPLEX(xp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: moments_save
  COMPLEX(xp), DIMENSION(:,:,:),       ALLOCATABLE :: phi_save, psi_save, dBpar_save
  ! Chunk of blocks currently loaded in the compact Krylov representation
  INTEGER :: nblk_c = 0, nmode_c = 0, ncomp_c = 0
  INTEGER, DIMENSION(:), ALLOCATABLE :: mode_iky, mode_ikx, blk_ms, blk_me, blk_bs, blk_be

  ABSTRACT INTERFACE
    SUBROUTINE matvec_iface(v, w)
      IMPORT :: xp
      COMPLEX(xp), DIMENSION(:), INTENT(IN)  :: v
      COMPLEX(xp), DIMENSION(:), INTENT(OUT) :: w
    END SUBROUTINE matvec_iface
  END INTERFACE

  PUBLIC :: adaptive_dt_readinputs, adaptive_dt_outputinputs, adaptive_dt_finalize_output
  PUBLIC :: adaptive_dt_init, adaptive_dt_update
  ! Building blocks, public for the unit tests
  PUBLIC :: matvec_iface, setup_stability, stability_function, h_stab, h_of_omega
  PUBLIC :: hessenberg_eig, ritz_residuals, hermite_ladder_radius, arnoldi_blocks
  PUBLIC :: next_dt, compute_omega_nl, estimate_linear_limit, build_h_table

CONTAINS

  !=============================================================================
  !                               Input / output
  !=============================================================================
  SUBROUTINE adaptive_dt_readinputs
    ! Read the ADAPTIVE_DT namelist. The group is optional and may sit anywhere
    ! in the file (REWIND + iostat, like LINEAR_GROWTH).
    USE basic, ONLY: lu_in
    IMPLICIT NONE
    INTEGER :: iost
    NAMELIST /ADAPTIVE_DT/ use_adaptive_dt, cfl_frac, dt_min, dt_max, dt_growth, dt_hyst, &
                           estimate_dt_lin, dt_lin, dt_lin_safety, krylov_dim, krylov_mem_levels, &
                           n_omega_tab, dt_verbose
    REWIND(lu_in)
    READ(lu_in, nml=ADAPTIVE_DT, iostat=iost)
    IF (iost .GT. 0) ERROR STOP '>> ERROR << could not parse the ADAPTIVE_DT namelist'
    IF (.NOT. use_adaptive_dt) RETURN
    IF (cfl_frac .LE. 0._xp)          ERROR STOP '>> ERROR << ADAPTIVE_DT: cfl_frac must be positive'
    IF (dt_min .LE. 0._xp)            ERROR STOP '>> ERROR << ADAPTIVE_DT: dt_min must be positive'
    IF (dt_max .LT. dt_min)           ERROR STOP '>> ERROR << ADAPTIVE_DT: need dt_max >= dt_min'
    IF (dt_growth .LE. 1._xp)         ERROR STOP '>> ERROR << ADAPTIVE_DT: dt_growth must exceed 1'
    IF (dt_hyst .LT. 0._xp)           ERROR STOP '>> ERROR << ADAPTIVE_DT: dt_hyst must be >= 0'
    IF (krylov_dim .LT. 2 .OR. krylov_dim .GT. 128) &
                                      ERROR STOP '>> ERROR << ADAPTIVE_DT: krylov_dim must be in [2,128]'
    IF (krylov_mem_levels .LE. 0._xp) ERROR STOP '>> ERROR << ADAPTIVE_DT: krylov_mem_levels must be positive'
    IF (n_omega_tab .LT. 4)           ERROR STOP '>> ERROR << ADAPTIVE_DT: n_omega_tab must be >= 4'
    IF (dt_lin_safety .LE. 0._xp .OR. dt_lin_safety .GT. 1._xp) &
                                      ERROR STOP '>> ERROR << ADAPTIVE_DT: dt_lin_safety must be in (0,1]'
    IF ((.NOT. estimate_dt_lin) .AND. (dt_lin .LE. 0._xp)) &
      ERROR STOP '>> ERROR << ADAPTIVE_DT: give a positive dt_lin when estimate_dt_lin is off'
  END SUBROUTINE adaptive_dt_readinputs

  SUBROUTINE adaptive_dt_outputinputs(fid)
    ! Write the parameters and the result of the linear estimate to the output
    ! file. A group (not a 0d dataset) because it also holds the dt_lin map.
    USE futils,   ONLY: creatg, attach, putarr
    USE grid,     ONLY: local_nky, local_nkx, local_nz, total_nky, total_nkx, total_nz
    USE parallel, ONLY: gather_xyz_real, num_procs
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fid
    CHARACTER(len=256)  :: str
    REAL(xp), DIMENSION(:,:,:), ALLOCATABLE :: map_loc, map_full
    INTEGER :: iz
    WRITE(str,'(a)') '/data/input/adaptive_dt'
    CALL creatg(fid, TRIM(str), 'Adaptive time step input and linear stability estimate')
    CALL attach(fid, TRIM(str), "use_adaptive_dt",   use_adaptive_dt)
    CALL attach(fid, TRIM(str), "cfl_frac",          cfl_frac)
    CALL attach(fid, TRIM(str), "dt_min",            dt_min)
    CALL attach(fid, TRIM(str), "dt_max",            dt_max)
    CALL attach(fid, TRIM(str), "dt_growth",         dt_growth)
    CALL attach(fid, TRIM(str), "dt_hyst",           dt_hyst)
    CALL attach(fid, TRIM(str), "estimate_dt_lin",   estimate_dt_lin)
    CALL attach(fid, TRIM(str), "dt_lin",            dt_lin)
    CALL attach(fid, TRIM(str), "dt_lin_safety",     dt_lin_safety)
    CALL attach(fid, TRIM(str), "krylov_dim",        krylov_dim)
    CALL attach(fid, TRIM(str), "krylov_mem_levels", krylov_mem_levels)
    CALL attach(fid, TRIM(str), "n_omega_tab",       n_omega_tab)
    CALL attach(fid, TRIM(str), "S_real",            S_real)
    CALL attach(fid, TRIM(str), "S_imag",            S_imag)
    CALL attach(fid, TRIM(str), "row_blocks",        row_blocks)
    CALL attach(fid, TRIM(str), "lim_kx",            lim_kx)
    CALL attach(fid, TRIM(str), "lim_ky",            lim_ky)
    CALL attach(fid, TRIM(str), "lim_eig_re",        REAL(lim_eig,xp))
    CALL attach(fid, TRIM(str), "lim_eig_im",        AIMAG(lim_eig))
    CALL attach(fid, TRIM(str), "lim_resid",         lim_resid)
    CALL attach(fid, TRIM(str), "n_rhs_evals",       n_rhs_evals)
    CALL attach(fid, TRIM(str), "arnoldi_time",      arnoldi_time)
    CALL attach(fid, TRIM(str), "dt_initial",        dt_used)
    IF (ALLOCATED(dt_lin_map)) THEN
      ! gather_xyz_real works on (ky,kx,z) fields: replicate the map along z
      ALLOCATE(map_loc(local_nky,local_nkx,local_nz), map_full(total_nky,total_nkx,total_nz))
      DO iz = 1,local_nz
        map_loc(:,:,iz) = dt_lin_map
      ENDDO
      map_full = 0._xp
      IF (num_procs .EQ. 1) THEN
        map_full = map_loc
      ELSE
        CALL gather_xyz_real(map_loc, map_full, local_nky, total_nky, total_nkx, local_nz, total_nz)
      ENDIF
      CALL putarr(fid, TRIM(str)//'/dt_lin_kykx', map_full(:,:,1), 'linear stability limit per (ky,kx) mode', ionode=0)
      DEALLOCATE(map_loc, map_full)
    ENDIF
  END SUBROUTINE adaptive_dt_outputinputs

  SUBROUTINE adaptive_dt_finalize_output(fid)
    ! Counters of the whole run, attached at the last diagnostic.
    USE futils, ONLY: attach
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fid
    CALL attach(fid, '/data/input/adaptive_dt', "n_dt_changes",     n_dt_changes)
    CALL attach(fid, '/data/input/adaptive_dt', "n_cfl_violations", n_cfl_violations)
    CALL attach(fid, '/data/input/adaptive_dt', "n_cfl_warnings",   n_cfl_warnings)
    CALL attach(fid, '/data/input/adaptive_dt', "max_cfl_excess",   max_cfl_excess)
  END SUBROUTINE adaptive_dt_finalize_output

  !=============================================================================
  !                           Driver: init and update
  !=============================================================================
  SUBROUTINE adaptive_dt_init
    ! Called once after the initial state is set and before the first
    ! diagnostic: measures the limits and sets the first time step.
    USE basic,            ONLY: speak, time, tmax, change_dt, dt
    USE model,            ONLY: EM
    USE grid,             ONLY: pmax
    USE species,          ONLY: sqrt_tau_o_sigma
    USE time_integration, ONLY: numerical_scheme
    IMPLICIT NONE
    CHARACTER(len=256) :: str
    REAL(xp) :: dt0
    IF (.NOT. use_adaptive_dt) RETURN
    CALL speak('Adaptive time step [',1)
    CALL setup_stability
    WRITE(str,'(A,A,A,F7.4,A,F7.4)') '-scheme ', TRIM(numerical_scheme), &
      ': stability interval S_real = ', S_real, ', S_imag = ', S_imag
    CALL speak(TRIM(str),0)
    IF (S_imag .LT. 1.0e-3_xp) ERROR STOP &
      '>> ERROR << ADAPTIVE_DT needs a scheme with a stability interval on the imaginary axis (RK3, SSP_RK3, ODE23, RK4, DOPRI5)'
    flutter_factor = 0._xp
    IF (EM) flutter_factor = MAXVAL(sqrt_tau_o_sigma)*hermite_ladder_radius(pmax)
    IF (estimate_dt_lin) THEN
      CALL estimate_linear_limit
      WRITE(str,'(A,ES11.4,A,F8.3,A,F8.3,A,ES11.3,A,ES11.3,A,ES8.1,A)') '-linear limit dt_lin = ', dt_lin, &
        '  (ky = ', lim_ky, ', kx = ', lim_kx, ', eig = ', REAL(lim_eig,xp), ' + i ', AIMAG(lim_eig), &
        ', resid = ', lim_resid, ')'
      CALL speak(TRIM(str),0)
      WRITE(str,'(A,I0,A,I0,A,F8.2,A)') '-Arnoldi: krylov_dim = ', krylov_dim, ', ', n_rhs_evals, &
        ' RHS evaluations, ', arnoldi_time, ' s'
      CALL speak(TRIM(str),0)
    ELSE
      WRITE(str,'(A,ES11.4)') '-linear limit dt_lin (input) = ', dt_lin
      CALL speak(TRIM(str),0)
    ENDIF
    CALL build_h_table
    CALL compute_omega_nl
    dt_est = h_of_omega(omega_nl)
    dt0 = next_dt(HUGE_H, cfl_frac*dt_est, time, tmax) ! HUGE: no growth limit on the first value
    IF (cfl_frac*dt_est .LT. dt_min) ERROR STOP '>> ERROR << ADAPTIVE_DT: the initial stable time step is below dt_min'
    CALL change_dt(dt0)
    dt_used = dt0
    WRITE(str,'(A,ES11.4,A,ES11.4,A,ES11.4,A)') '-initial dt = ', dt, '  (omega_nl = ', omega_nl, &
      ', stability limit = ', dt_est, ')'
    CALL speak(TRIM(str),0)
    CALL speak('] adaptive time step ready',1)
  END SUBROUTINE adaptive_dt_init

  SUBROUTINE adaptive_dt_update
    ! Called at the end of every iteration of the main loop (after the
    ! diagnostics, so that the outputs carry the dt that produced the state).
    ! Measures the CFL of the new state and sets dt for the next step.
    USE basic,    ONLY: dt, time, tmax, change_dt, speak
    USE parallel, ONLY: my_id
    USE mpi
    IMPLICIT NONE
    REAL(xp) :: dt_target, dt_new
    INTEGER  :: ierr
    CHARACTER(len=256) :: str
    IF (.NOT. use_adaptive_dt) RETURN
    dt_used = dt
    CALL compute_omega_nl
    dt_est = h_of_omega(omega_nl)
    ! The step just taken is checked against the limit the new state implies.
    ! The velocity is sampled at the start of a step, so exceedances of a few
    ! percent are inherent and harmless (the amplification just past the
    ! boundary is marginal); they are counted, while the warning is kept for
    ! exceedances that suggest a cfl_frac above what the run tolerates.
    IF (dt_used .GT. dt_est*(1._xp + 0.01_xp)) THEN
      n_cfl_violations = n_cfl_violations + 1
      max_cfl_excess   = MAX(max_cfl_excess, dt_used/dt_est - 1._xp)
      IF (dt_used .GT. dt_est*1.1_xp) THEN
        n_cfl_warnings = n_cfl_warnings + 1
        IF (n_cfl_warnings .EQ. 1 .OR. MOD(n_cfl_warnings,100) .EQ. 0) THEN
          WRITE(str,'(A,ES10.3,A,ES10.3,A,I0,A)') 'ADAPTIVE_DT: dt = ', dt_used, &
            ' exceeded the CFL limit ', dt_est, ' by more than 10% (', n_cfl_warnings, &
            ' times); consider lowering cfl_frac'
          CALL speak(TRIM(str),0)
        ENDIF
      ENDIF
    ENDIF
    dt_target = cfl_frac*dt_est
    dt_new    = next_dt(dt, dt_target, time, tmax)
    ! every rank computed the same value from reduced data; make it bitwise so
    CALL MPI_BCAST(dt_new, 1, mpi_xp_r, 0, MPI_COMM_WORLD, ierr)
    IF (dt_target .LT. dt_min) THEN
      adt_stop = .TRUE. ! tesend ends the run cleanly on this flag
      WRITE(str,'(A,ES10.3,A,ES10.3,A)') 'ADAPTIVE_DT: stable time step ', dt_target, &
        ' fell below dt_min = ', dt_min, ' -> stopping'
      CALL speak(TRIM(str),-99)
      RETURN
    ENDIF
    IF (dt_new .NE. dt) THEN
      n_dt_changes = n_dt_changes + 1
      IF (dt_verbose .AND. my_id .EQ. 0) THEN
        WRITE(str,'(A,F10.3,A,ES10.3,A,ES10.3,A,ES10.3,A)') 'ADAPTIVE_DT: t = ', time, &
          ' dt ', dt, ' -> ', dt_new, ' (omega_nl = ', omega_nl, ')'
        CALL speak(TRIM(str),0)
      ENDIF
      CALL change_dt(dt_new)
    ENDIF
  END SUBROUTINE adaptive_dt_update

  !> The step size controller: immediate decrease, growth limited to dt_growth
  !! per step and only once the target exceeds the current dt by dt_hyst,
  !! bounds, and an exact landing on tmax.
  PURE FUNCTION next_dt(dt_now, dt_target, time_now, tmax_) RESULT(dt_new)
    IMPLICIT NONE
    REAL(xp), INTENT(IN) :: dt_now, dt_target, time_now, tmax_
    REAL(xp) :: dt_new, remaining
    dt_new = dt_now
    IF (dt_target .LT. dt_now) THEN
      dt_new = dt_target
    ELSE IF (dt_target .GT. (1._xp + dt_hyst)*dt_now) THEN
      dt_new = MIN(dt_target, dt_growth*dt_now)
    ENDIF
    dt_new = MIN(MAX(dt_new, dt_min), dt_max)
    remaining = tmax_ - time_now
    IF ((remaining .GT. 0._xp) .AND. (dt_new .GT. remaining)) &
      dt_new = remaining*(1._xp + 4._xp*EPSILON(1._xp)) ! lands at or just past tmax
  END FUNCTION next_dt

  !=============================================================================
  !                   Stability region of the Butcher tableau
  !=============================================================================
  !> Stability polynomial R(z) of the scheme, evaluated by running the scheme
  !! itself on y' = lambda y with z = h lambda: stage values
  !! g_i = 1 + z sum_j A_ij g_j and R = 1 + z sum_i b_i g_i. Generic for every
  !! explicit tableau of time_integration_mod and guaranteed to match it.
  FUNCTION stability_function(z) RESULT(R)
    USE time_integration, ONLY: A_E, b_E, ntimelevel
    IMPLICIT NONE
    COMPLEX(xp), INTENT(IN) :: z
    COMPLEX(xp) :: R
    COMPLEX(xp), DIMENSION(ntimelevel) :: g
    INTEGER :: i, j
    DO i = 1,ntimelevel
      g(i) = (1._xp, 0._xp)
      DO j = 1,i-1
        g(i) = g(i) + z*A_E(i,j)*g(j)
      ENDDO
    ENDDO
    R = (1._xp, 0._xp)
    DO i = 1,ntimelevel
      R = R + z*b_E(i,1)*g(i)
    ENDDO
  END FUNCTION stability_function

  LOGICAL FUNCTION is_stable(z)
    IMPLICIT NONE
    COMPLEX(xp), INTENT(IN) :: z
    is_stable = ABS(stability_function(z)) .LE. 1._xp + 10._xp*EPSILON(1._xp)
  END FUNCTION is_stable

  !> Extent of the stability region along the direction u (|u| = 1): the
  !! largest h such that h' u is stable for all h' in [0,h]. Slow but robust
  !! upward scan, used only to characterize the scheme.
  FUNCTION region_extent(u) RESULT(h)
    IMPLICIT NONE
    COMPLEX(xp), INTENT(IN) :: u
    REAL(xp) :: h, h_lo, h_hi
    INTEGER  :: it
    h_lo = 0._xp
    h    = 1.0e-3_xp
    DO WHILE (h .LT. 1.0e3_xp)
      IF (.NOT. is_stable(h*u)) EXIT
      h_lo = h
      h    = h*1.05_xp
    ENDDO
    IF (h .GE. 1.0e3_xp) THEN
      h = HUGE_H
      RETURN
    ENDIF
    h_hi = h
    DO it = 1,50
      h = 0.5_xp*(h_lo + h_hi)
      IF (is_stable(h*u)) THEN
        h_lo = h
      ELSE
        h_hi = h
      ENDIF
    ENDDO
    h = h_lo
  END FUNCTION region_extent

  !> Characterize the scheme: S_real, S_imag and the radius R_max of the
  !! smallest disk centered at 0 containing the stability region.
  SUBROUTINE setup_stability
    USE prec_const, ONLY: PI
    IMPLICIT NONE
    INTEGER :: k
    REAL(xp) :: theta, r
    S_real = region_extent((-1._xp, 0._xp))
    S_imag = region_extent(( 0._xp, 1._xp))
    R_max  = 0._xp
    DO k = 0,143
      theta = PI*(0.5_xp + REAL(k,xp)/144._xp) ! left half plane and the imaginary axis
      r = region_extent(CMPLX(COS(theta), SIN(theta), xp))
      IF (r .LT. HUGE_H) R_max = MAX(R_max, r)
    ENDDO
    R_max = MAX(R_max, S_real, S_imag)
  END SUBROUTINE setup_stability

  !> Largest stable step for the eigenvalue lambda: max h with |R(h' lambda)|
  !! <= 1 for all h' in [0,h]. Eigenvalues with a positive real part (physical
  !! instabilities, limited by accuracy and not by stability) are projected on
  !! the imaginary axis first. Requires setup_stability.
  FUNCTION h_stab(lambda) RESULT(h)
    IMPLICIT NONE
    COMPLEX(xp), INTENT(IN) :: lambda
    REAL(xp) :: h, h_lo, h_hi, alam
    COMPLEX(xp) :: lam
    INTEGER :: it
    lam  = CMPLX(MIN(REAL(lambda,xp), 0._xp), AIMAG(lambda), xp)
    alam = ABS(lam)
    IF (alam .LE. TINY(1._xp)) THEN
      h = HUGE_H
      RETURN
    ENDIF
    ! h_hi is unstable by containment in the disk of radius R_max
    h_hi = 1.01_xp*R_max/alam
    h_lo = 0._xp
    h    = h_hi
    DO it = 1,60
      h = 0.75_xp*h
      IF (is_stable(h*lam)) THEN
        h_lo = h
        EXIT
      ENDIF
      h_hi = h
    ENDDO
    IF (h_lo .EQ. 0._xp) THEN ! never stable along this direction (e.g. RK2 on the imaginary axis)
      h = 0._xp
      RETURN
    ENDIF
    DO it = 1,40
      h = 0.5_xp*(h_lo + h_hi)
      IF (is_stable(h*lam)) THEN
        h_lo = h
      ELSE
        h_hi = h
      ENDIF
    ENDDO
    h = h_lo
  END FUNCTION h_stab

  !=============================================================================
  !                 Small dense linear algebra (no LAPACK needed)
  !=============================================================================
  !> Eigenvalues of an n x n complex upper Hessenberg matrix by the shifted QR
  !! algorithm with Wilkinson shifts and deflation (Givens rotations). info is
  !! the number of eigenvalues that did not converge (0 normally).
  SUBROUTINE hessenberg_eig(n, Hin, eig, info)
    IMPLICIT NONE
    INTEGER,                     INTENT(IN)  :: n
    COMPLEX(xp), DIMENSION(n,n), INTENT(IN)  :: Hin
    COMPLEX(xp), DIMENSION(n),   INTENT(OUT) :: eig
    INTEGER,                     INTENT(OUT) :: info
    COMPLEX(xp), DIMENSION(n,n) :: B
    COMPLEX(xp), DIMENSION(n)   :: cs, sn
    COMPLEX(xp) :: e11, e12, e21, e22, mu, mu1, mu2, disc, t1, t2, cc, ss
    REAL(xp)    :: r, eps, hnorm
    INTEGER     :: hi, lo, k, i, j, iter, total_iter
    eps   = EPSILON(1._xp)
    B     = Hin
    hnorm = SQRT(SUM(ABS(B)**2))
    info  = 0
    eig   = (0._xp, 0._xp)
    hi = n; iter = 0; total_iter = 0
    DO WHILE (hi .GE. 1)
      ! active block is lo..hi, lo found from a negligible subdiagonal
      lo = hi
      DO WHILE (lo .GT. 1)
        IF (ABS(B(lo,lo-1)) .LE. eps*(ABS(B(lo-1,lo-1)) + ABS(B(lo,lo))) .OR. &
            ABS(B(lo,lo-1)) .LE. 10._xp*eps*hnorm) EXIT
        lo = lo - 1
      ENDDO
      IF (lo .EQ. hi) THEN
        eig(hi) = B(hi,hi)
        hi = hi - 1
        iter = 0
        CYCLE
      ENDIF
      iter = iter + 1; total_iter = total_iter + 1
      IF (total_iter .GT. 200*n) THEN
        info = hi
        DO i = 1,hi
          eig(i) = B(i,i)
        ENDDO
        EXIT
      ENDIF
      ! shift: eigenvalue of the trailing 2x2 closest to B(hi,hi), or an
      ! exceptional one every so often to break cycles
      e11 = B(hi-1,hi-1); e12 = B(hi-1,hi); e21 = B(hi,hi-1); e22 = B(hi,hi)
      IF (MOD(iter,11) .EQ. 0) THEN
        mu = e22 + ABS(e21)*(1._xp + 0.5_xp*imagu)
      ELSE
        disc = SQRT(((e11-e22)*0.5_xp)**2 + e12*e21)
        mu1  = (e11+e22)*0.5_xp + disc
        mu2  = (e11+e22)*0.5_xp - disc
        IF (ABS(mu1-e22) .LT. ABS(mu2-e22)) THEN
          mu = mu1
        ELSE
          mu = mu2
        ENDIF
      ENDIF
      ! QR step on the active block: (B - mu I) = QR, B <- RQ + mu I
      DO k = lo,hi
        B(k,k) = B(k,k) - mu
      ENDDO
      DO k = lo,hi-1
        e11 = B(k,k); e21 = B(k+1,k)
        r = SQRT(ABS(e11)**2 + ABS(e21)**2)
        IF (r .EQ. 0._xp) THEN
          cc = (1._xp, 0._xp); ss = (0._xp, 0._xp)
        ELSE
          cc = e11/r; ss = e21/r
        ENDIF
        cs(k) = cc; sn(k) = ss
        DO j = k,hi
          t1 = B(k,j); t2 = B(k+1,j)
          B(k,j)   =  CONJG(cc)*t1 + CONJG(ss)*t2
          B(k+1,j) = -ss*t1 + cc*t2
        ENDDO
      ENDDO
      DO k = lo,hi-1
        cc = cs(k); ss = sn(k)
        DO i = lo,MIN(k+1,hi)
          t1 = B(i,k); t2 = B(i,k+1)
          B(i,k)   =  cc*t1 + ss*t2
          B(i,k+1) = -CONJG(ss)*t1 + CONJG(cc)*t2
        ENDDO
      ENDDO
      DO k = lo,hi
        B(k,k) = B(k,k) + mu
      ENDDO
    ENDDO
  END SUBROUTINE hessenberg_eig

  !> LU factorization with partial pivoting of a dense complex matrix. Exact
  !! zero pivots are replaced by a tiny value (the matrices factored here are
  !! near singular by construction, see ritz_residuals).
  SUBROUTINE lu_factor(n, A, piv, pivtiny)
    IMPLICIT NONE
    INTEGER,                     INTENT(IN)    :: n
    COMPLEX(xp), DIMENSION(n,n), INTENT(INOUT) :: A
    INTEGER,     DIMENSION(n),   INTENT(OUT)   :: piv
    REAL(xp),                    INTENT(IN)    :: pivtiny
    INTEGER :: k, i, p
    COMPLEX(xp), DIMENSION(n) :: tmp
    DO k = 1,n
      p = k
      DO i = k+1,n
        IF (ABS(A(i,k)) .GT. ABS(A(p,k))) p = i
      ENDDO
      piv(k) = p
      IF (p .NE. k) THEN
        tmp = A(k,:); A(k,:) = A(p,:); A(p,:) = tmp
      ENDIF
      IF (ABS(A(k,k)) .LE. pivtiny) A(k,k) = CMPLX(pivtiny, 0._xp, xp)
      DO i = k+1,n
        A(i,k) = A(i,k)/A(k,k)
        A(i,k+1:n) = A(i,k+1:n) - A(i,k)*A(k,k+1:n)
      ENDDO
    ENDDO
  END SUBROUTINE lu_factor

  SUBROUTINE lu_solve(n, A, piv, b)
    IMPLICIT NONE
    INTEGER,                     INTENT(IN)    :: n
    COMPLEX(xp), DIMENSION(n,n), INTENT(IN)    :: A
    INTEGER,     DIMENSION(n),   INTENT(IN)    :: piv
    COMPLEX(xp), DIMENSION(n),   INTENT(INOUT) :: b
    INTEGER :: k, i
    COMPLEX(xp) :: t
    DO k = 1,n
      IF (piv(k) .NE. k) THEN
        t = b(k); b(k) = b(piv(k)); b(piv(k)) = t
      ENDIF
      DO i = k+1,n
        b(i) = b(i) - A(i,k)*b(k)
      ENDDO
    ENDDO
    DO k = n,1,-1
      b(k) = b(k)/A(k,k)
      DO i = 1,k-1
        b(i) = b(i) - A(i,k)*b(k)
      ENDDO
    ENDDO
  END SUBROUTINE lu_solve

  !> Residual norms ||L x - lambda x|| = |h_{m+1,m}| |y_m| of the Ritz pairs of
  !! an m-step Arnoldi factorization, with the eigenvectors y of the Hessenberg
  !! matrix obtained by two steps of inverse iteration.
  SUBROUTINE ritz_residuals(m, H, hsub, eig, resid)
    IMPLICIT NONE
    INTEGER,                     INTENT(IN)  :: m
    COMPLEX(xp), DIMENSION(m,m), INTENT(IN)  :: H
    REAL(xp),                    INTENT(IN)  :: hsub
    COMPLEX(xp), DIMENSION(m),   INTENT(IN)  :: eig
    REAL(xp),    DIMENSION(m),   INTENT(OUT) :: resid
    COMPLEX(xp), DIMENSION(m,m) :: A
    COMPLEX(xp), DIMENSION(m)   :: y
    INTEGER,     DIMENSION(m)   :: piv
    REAL(xp) :: hnorm, pivtiny, ynorm
    INTEGER  :: k, i, it
    hnorm = SQRT(SUM(ABS(H)**2))
    IF (hnorm .LE. TINY(1._xp)) THEN
      ! a zero block (nothing acts on these modes): exact, with no residual
      resid = 0._xp
      RETURN
    ENDIF
    pivtiny = 10._xp*EPSILON(1._xp)*hnorm
    DO k = 1,m
      A = H
      DO i = 1,m
        A(i,i) = A(i,i) - eig(k)
      ENDDO
      CALL lu_factor(m, A, piv, pivtiny)
      y = (1._xp, 0._xp)
      DO it = 1,2
        CALL lu_solve(m, A, piv, y)
        ynorm = SQRT(SUM(ABS(y)**2))
        IF (ynorm .GT. 0._xp) y = y/ynorm
      ENDDO
      resid(k) = hsub*ABS(y(m))
    ENDDO
  END SUBROUTINE ritz_residuals

  !> Spectral radius of the Hermite ladder: the symmetric tridiagonal matrix
  !! with off diagonals sqrt(p+1), p = 0..P-1, whose eigenvalues are the zeros
  !! of He_{P+1}. It bounds the parallel velocity carried by moments up to P,
  !! i.e. the magnetic flutter advection speed of the sqrt(p+1)/sqrt(p) bracket.
  FUNCTION hermite_ladder_radius(P) RESULT(rad)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: P
    REAL(xp) :: rad
    COMPLEX(xp), DIMENSION(P+1,P+1) :: J
    COMPLEX(xp), DIMENSION(P+1)     :: eig
    INTEGER :: ip, info
    IF (P .LT. 1) THEN
      rad = 0._xp
      RETURN
    ENDIF
    J = (0._xp, 0._xp)
    DO ip = 1,P
      J(ip,ip+1) = CMPLX(SQRT(REAL(ip,xp)), 0._xp, xp)
      J(ip+1,ip) = J(ip,ip+1)
    ENDDO
    CALL hessenberg_eig(P+1, J, eig, info)
    rad = MAXVAL(ABS(eig))
  END FUNCTION hermite_ladder_radius

  !=============================================================================
  !                       Block Arnoldi on compact vectors
  !=============================================================================
  !> Deterministic pseudo random start vector (no RANDOM_NUMBER, so the state
  !! of the generator used by the initial conditions is untouched).
  SUBROUTINE lcg_fill(v, seed)
    ! Park-Miller minimal standard generator in 64-bit integers (no overflow)
    IMPLICIT NONE
    COMPLEX(xp), DIMENSION(:), INTENT(OUT) :: v
    INTEGER, INTENT(IN) :: seed
    INTEGER(KIND=8), PARAMETER :: a = 48271_8, mm = 2147483647_8
    INTEGER(KIND=8) :: s
    INTEGER :: i
    REAL(xp) :: u1, u2
    s = MODULO(INT(seed,8), mm)
    IF (s .EQ. 0_8) s = 1_8
    DO i = 1,SIZE(v)
      s  = MOD(a*s, mm)
      u1 = REAL(s,xp)/REAL(mm,xp)
      s  = MOD(a*s, mm)
      u2 = REAL(s,xp)/REAL(mm,xp)
      v(i) = CMPLX(u1 - 0.5_xp, u2 - 0.5_xp, xp)
    ENDDO
  END SUBROUTINE lcg_fill

  !> Arnoldi iterations run simultaneously on nblk independent blocks of a
  !! block diagonal operator. The compact vectors of length ncomp hold the
  !! blocks one after the other (block ib occupies bs(ib):be(ib)); every inner
  !! product is taken per block and reduced over comm (the ranks sharing the
  !! blocks). All ranks perform exactly m matvecs; a block whose Krylov space
  !! becomes invariant (breakdown) is frozen, and m_eff holds the dimension it
  !! reached. H(1:m_eff+1,1:m_eff,ib) is the Hessenberg matrix of block ib.
  SUBROUTINE arnoldi_blocks(matvec, ncomp, nblk, bs, be, m, comm, seed, H, m_eff)
    USE mpi
    IMPLICIT NONE
    PROCEDURE(matvec_iface) :: matvec
    INTEGER, INTENT(IN) :: ncomp, nblk, m, comm, seed
    INTEGER, DIMENSION(nblk), INTENT(IN) :: bs, be
    COMPLEX(xp), DIMENSION(m+1,m,nblk), INTENT(OUT) :: H
    INTEGER, DIMENSION(nblk), INTENT(OUT) :: m_eff
    COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE :: V, dots
    COMPLEX(xp), DIMENSION(:),   ALLOCATABLE :: w
    REAL(xp),    DIMENSION(:),   ALLOCATABLE :: nrm
    LOGICAL,     DIMENSION(:),   ALLOCATABLE :: active
    INTEGER  :: i, j, ib, ierr, pass
    REAL(xp) :: colnorm, breakdown_tol
    ALLOCATE(V(ncomp,m+1), w(ncomp), dots(m,nblk), nrm(nblk), active(nblk))
    H = (0._xp, 0._xp); m_eff = 0
    breakdown_tol = 100._xp*EPSILON(1._xp)
    ! start vector, normalized per block
    CALL lcg_fill(V(:,1), seed)
    CALL block_norms(V(:,1), nrm)
    DO ib = 1,nblk
      active(ib) = nrm(ib) .GT. 0._xp
      IF (active(ib)) THEN
        V(bs(ib):be(ib),1) = V(bs(ib):be(ib),1)/SQRT(nrm(ib))
      ELSE
        V(bs(ib):be(ib),1) = (0._xp, 0._xp)
      ENDIF
    ENDDO
    DO j = 1,m
      CALL matvec(V(:,j), w)
      ! modified Gram-Schmidt, twice (the second pass restores orthogonality
      ! lost to round off, which matters in the single precision build)
      DO pass = 1,2
        dots = (0._xp, 0._xp)
        DO ib = 1,nblk
          IF (.NOT. active(ib)) CYCLE
          DO i = 1,j
            dots(i,ib) = SUM(CONJG(V(bs(ib):be(ib),i))*w(bs(ib):be(ib)))
          ENDDO
        ENDDO
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, dots, m*nblk, mpi_xp_c, MPI_SUM, comm, ierr)
        DO ib = 1,nblk
          IF (.NOT. active(ib)) CYCLE
          DO i = 1,j
            w(bs(ib):be(ib)) = w(bs(ib):be(ib)) - dots(i,ib)*V(bs(ib):be(ib),i)
            H(i,j,ib) = H(i,j,ib) + dots(i,ib)
          ENDDO
        ENDDO
      ENDDO
      CALL block_norms(w, nrm)
      DO ib = 1,nblk
        V(bs(ib):be(ib),j+1) = (0._xp, 0._xp)
        IF (.NOT. active(ib)) CYCLE
        H(j+1,j,ib) = CMPLX(SQRT(nrm(ib)), 0._xp, xp)
        colnorm = SQRT(SUM(ABS(H(1:j+1,j,ib))**2))
        m_eff(ib) = j
        IF (REAL(H(j+1,j,ib),xp) .LE. breakdown_tol*colnorm) THEN
          ! invariant subspace found: the Ritz values of this block are exact
          H(j+1,j,ib) = (0._xp, 0._xp)
          active(ib)  = .FALSE.
        ELSE
          V(bs(ib):be(ib),j+1) = w(bs(ib):be(ib))/REAL(H(j+1,j,ib),xp)
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(V, w, dots, nrm, active)
  CONTAINS
    SUBROUTINE block_norms(x, nrm_)
      COMPLEX(xp), DIMENSION(:), INTENT(IN)  :: x
      REAL(xp),    DIMENSION(:), INTENT(OUT) :: nrm_
      INTEGER :: ib_, ierr_
      DO ib_ = 1,nblk
        nrm_(ib_) = SUM(ABS(x(bs(ib_):be(ib_)))**2)
      ENDDO
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, nrm_, nblk, mpi_xp_r, MPI_SUM, comm, ierr_)
    END SUBROUTINE block_norms
  END SUBROUTINE arnoldi_blocks

  !=============================================================================
  !        Linear operator of the run, matrix free, on compact vectors
  !=============================================================================
  !> Define the blocks of chunk ic (0 <= ic < nchunks) on this rank and the
  !! compact layout of their evolved entries.
  SUBROUTINE setup_chunk(ic, nchunks)
    USE grid,    ONLY: local_nky, local_nkx, local_nz, local_na, local_np, local_nj, &
                       ngp_o2, ngj_o2, AA_x, AA_y, local_nky_offset, total_nkx
    USE closure, ONLY: evolve_mom
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ic, nchunks
    INTEGER :: iky, ikx, nev, ib, im, gid
    nev = local_na*COUNT(evolve_mom(1+ngp_o2:local_np+ngp_o2, 1+ngj_o2:local_nj+ngj_o2))
    ! count first
    nblk_c = 0; nmode_c = 0
    DO iky = 1,local_nky
      IF (AA_y(iky) .LE. 0._xp) CYCLE
      IF (row_blocks) THEN
        gid = iky + local_nky_offset
        IF (MOD(gid, nchunks) .NE. ic) CYCLE
        IF (COUNT(AA_x .GT. 0._xp) .EQ. 0) CYCLE
        nblk_c  = nblk_c + 1
        nmode_c = nmode_c + COUNT(AA_x .GT. 0._xp)
      ELSE
        DO ikx = 1,local_nkx
          IF (AA_x(ikx) .LE. 0._xp) CYCLE
          gid = (iky + local_nky_offset - 1)*total_nkx + ikx
          IF (MOD(gid, nchunks) .NE. ic) CYCLE
          nblk_c  = nblk_c + 1
          nmode_c = nmode_c + 1
        ENDDO
      ENDIF
    ENDDO
    IF (ALLOCATED(mode_iky)) DEALLOCATE(mode_iky, mode_ikx, blk_ms, blk_me, blk_bs, blk_be)
    ALLOCATE(mode_iky(MAX(nmode_c,1)), mode_ikx(MAX(nmode_c,1)))
    ALLOCATE(blk_ms(MAX(nblk_c,1)), blk_me(MAX(nblk_c,1)), blk_bs(MAX(nblk_c,1)), blk_be(MAX(nblk_c,1)))
    ! then fill
    ib = 0; im = 0; ncomp_c = 0
    DO iky = 1,local_nky
      IF (AA_y(iky) .LE. 0._xp) CYCLE
      IF (row_blocks) THEN
        gid = iky + local_nky_offset
        IF (MOD(gid, nchunks) .NE. ic) CYCLE
        IF (COUNT(AA_x .GT. 0._xp) .EQ. 0) CYCLE
        ib = ib + 1
        blk_ms(ib) = im + 1
        DO ikx = 1,local_nkx
          IF (AA_x(ikx) .LE. 0._xp) CYCLE
          im = im + 1
          mode_iky(im) = iky; mode_ikx(im) = ikx
        ENDDO
        blk_me(ib) = im
        blk_bs(ib) = ncomp_c + 1
        ncomp_c    = ncomp_c + (blk_me(ib) - blk_ms(ib) + 1)*local_nz*nev
        blk_be(ib) = ncomp_c
      ELSE
        DO ikx = 1,local_nkx
          IF (AA_x(ikx) .LE. 0._xp) CYCLE
          gid = (iky + local_nky_offset - 1)*total_nkx + ikx
          IF (MOD(gid, nchunks) .NE. ic) CYCLE
          ib = ib + 1; im = im + 1
          mode_iky(im) = iky; mode_ikx(im) = ikx
          blk_ms(ib) = im; blk_me(ib) = im
          blk_bs(ib) = ncomp_c + 1
          ncomp_c    = ncomp_c + local_nz*nev
          blk_be(ib) = ncomp_c
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE setup_chunk

  !> Compact vector -> evolved interior entries of moments(...,1) of the
  !! loaded chunk (the rest of the level is kept at zero by the caller).
  SUBROUTINE scatter_to_moments(v)
    USE fields,  ONLY: moments
    USE grid,    ONLY: local_nz, local_na, local_np, local_nj, ngp_o2, ngj_o2, ngz_o2
    USE closure, ONLY: evolve_mom
    IMPLICIT NONE
    COMPLEX(xp), DIMENSION(:), INTENT(IN) :: v
    INTEGER :: ib, im, iky, ikx, iz, ij, ip, ia, izi, iji, ipi, ic
    ic = 0
    DO ib = 1,nblk_c
      DO im = blk_ms(ib),blk_me(ib)
        iky = mode_iky(im); ikx = mode_ikx(im)
        DO iz = 1,local_nz
          izi = iz + ngz_o2
          DO ij = 1,local_nj
            iji = ij + ngj_o2
            DO ip = 1,local_np
              ipi = ip + ngp_o2
              IF (.NOT. evolve_mom(ipi,iji)) CYCLE
              DO ia = 1,local_na
                ic = ic + 1
                moments(ia,ipi,iji,iky,ikx,izi,1) = v(ic)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE scatter_to_moments

  !> Evolved entries of moments_rhs(...,1) of the loaded chunk -> compact vector.
  SUBROUTINE gather_from_rhs(w)
    USE array,   ONLY: moments_rhs
    USE grid,    ONLY: local_nz, local_na, local_np, local_nj, ngp_o2, ngj_o2
    USE closure, ONLY: evolve_mom
    IMPLICIT NONE
    COMPLEX(xp), DIMENSION(:), INTENT(OUT) :: w
    INTEGER :: ib, im, iky, ikx, iz, ij, ip, ia, iji, ipi, ic
    ic = 0
    DO ib = 1,nblk_c
      DO im = blk_ms(ib),blk_me(ib)
        iky = mode_iky(im); ikx = mode_ikx(im)
        DO iz = 1,local_nz
          DO ij = 1,local_nj
            iji = ij + ngj_o2
            DO ip = 1,local_np
              ipi = ip + ngp_o2
              IF (.NOT. evolve_mom(ipi,iji)) CYCLE
              DO ia = 1,local_na
                ic = ic + 1
                w(ic) = moments_rhs(ia,ip,ij,iky,ikx,iz,1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE gather_from_rhs

  !> w = L v with L the linear operator of the run: the state goes through
  !! exactly the path of an RK stage (closure, ghosts, field solves, RHS) with
  !! the Poisson bracket left out.
  SUBROUTINE linear_rhs_matvec(v, w)
    USE stepper, ONLY: assemble_RHS, close_and_solve_fields
    IMPLICIT NONE
    COMPLEX(xp), DIMENSION(:), INTENT(IN)  :: v
    COMPLEX(xp), DIMENSION(:), INTENT(OUT) :: w
    CALL scatter_to_moments(v)
    CALL close_and_solve_fields
    CALL assemble_RHS(skip_nonlinear=.TRUE.)
    CALL gather_from_rhs(w)
    n_rhs_evals = n_rhs_evals + 1
  END SUBROUTINE linear_rhs_matvec

  !> Measure the linear stability limit: Arnoldi on every block, Ritz values
  !! pushed outward by their residual and by the safety margin, tested against
  !! the stability region. Fills dt_lin, dt_lin_map, the limiting block data
  !! and the list of extreme Ritz values used by the omega table. The state of
  !! the run is saved and restored around the pass.
  SUBROUTINE estimate_linear_limit
    USE fields,    ONLY: moments, phi, psi, dBpar
    USE grid,      ONLY: local_nky, local_nkx, ngz, AA_x, AA_y, kxarray, kyarray
    USE geometry,  ONLY: SHEARED, parallel_bc
    USE model,     ONLY: FORCE_SYMMETRY
    USE parallel,  ONLY: my_id, comm_pz, comm_ky
    USE nonlinear, ONLY: compute_Sapj
    USE mpi
    IMPLICIT NONE
    INTEGER  :: nblocks_global, nblocks_local, nchunks, ic, ib, im, k, me, info, nacc, ierr, m, owner
    INTEGER  :: nloc_cnt
    REAL(xp) :: t0, dt_lin_loc, hblk, hk, rel, best_rel
    REAL(xp), DIMENSION(4) :: lim_buf
    COMPLEX(xp), DIMENSION(:,:,:), ALLOCATABLE :: H
    COMPLEX(xp), DIMENSION(:,:),   ALLOCATABLE :: Hb
    COMPLEX(xp), DIMENSION(:),     ALLOCATABLE :: eig, lam
    REAL(xp),    DIMENSION(:),     ALLOCATABLE :: res, hlam
    INTEGER,     DIMENSION(:),     ALLOCATABLE :: m_eff
    LOGICAL,     DIMENSION(:),     ALLOCATABLE :: acc
    t0 = MPI_WTIME()
    n_rhs_evals = 0
    ! blocks: (kx,ky) modes, or whole ky rows when the parallel boundary
    ! condition chains kx modes (sheared geometry with z derivatives) or when
    ! the kx symmetry is enforced
    row_blocks = (SHEARED .AND. (ngz .GT. 0) .AND. (TRIM(parallel_bc) .NE. 'shearless') &
                  .AND. (TRIM(parallel_bc) .NE. 'disconnected')) .OR. FORCE_SYMMETRY
    IF (row_blocks) THEN
      nloc_cnt = 0
      IF (COUNT(AA_x .GT. 0._xp) .GT. 0) nloc_cnt = COUNT(AA_y .GT. 0._xp)
    ELSE
      nloc_cnt = COUNT(AA_y .GT. 0._xp)*COUNT(AA_x .GT. 0._xp)
    ENDIF
    nblocks_local = nloc_cnt
    CALL MPI_ALLREDUCE(nloc_cnt, nblocks_global, 1, MPI_INTEGER, MPI_SUM, comm_ky, ierr)
    IF (nblocks_global .LT. 1) ERROR STOP '>> ERROR << ADAPTIVE_DT: no evolved mode to estimate the linear limit on'
    m = krylov_dim
    nchunks = MAX(1, MIN(CEILING(REAL(m+1,xp)/krylov_mem_levels), nblocks_global))
    ! save the state
    ALLOCATE(moments_save(SIZE(moments,1),SIZE(moments,2),SIZE(moments,3),SIZE(moments,4),SIZE(moments,5),SIZE(moments,6)))
    ALLOCATE(phi_save(SIZE(phi,1),SIZE(phi,2),SIZE(phi,3)))
    ALLOCATE(psi_save(SIZE(psi,1),SIZE(psi,2),SIZE(psi,3)))
    ALLOCATE(dBpar_save(SIZE(dBpar,1),SIZE(dBpar,2),SIZE(dBpar,3)))
    moments_save = moments(:,:,:,:,:,:,1)
    phi_save = phi; psi_save = psi; dBpar_save = dBpar
    ! results
    IF (ALLOCATED(dt_lin_map)) DEALLOCATE(dt_lin_map)
    ALLOCATE(dt_lin_map(local_nky,local_nkx)); dt_lin_map = 0._xp
    IF (ALLOCATED(ritz_list)) DEALLOCATE(ritz_list)
    ALLOCATE(ritz_list(MAX(N_KEEP*nblocks_local,1))); n_ritz = 0
    ALLOCATE(eig(m), lam(m), res(m), hlam(m), acc(m))
    dt_lin_loc = HUGE_H
    lim_buf = 0._xp
    DO ic = 0,nchunks-1
      CALL setup_chunk(ic, nchunks)
      ALLOCATE(H(m+1,m,MAX(nblk_c,1)), m_eff(MAX(nblk_c,1)))
      moments(:,:,:,:,:,:,1) = (0._xp, 0._xp) ! everything outside the chunk is zero
      CALL arnoldi_blocks(linear_rhs_matvec, ncomp_c, nblk_c, blk_bs, blk_be, m, comm_pz, &
                          12345 + 7919*my_id + 104729*ic, H, m_eff)
      DO ib = 1,nblk_c
        me = m_eff(ib)
        IF (me .LT. 1) CYCLE
        ALLOCATE(Hb(me,me)); Hb = H(1:me,1:me,ib) ! contiguous copy for the explicit-shape routines
        CALL hessenberg_eig(me, Hb, eig(1:me), info)
        CALL ritz_residuals(me, Hb, ABS(H(me+1,me,ib)), eig(1:me), res(1:me))
        DEALLOCATE(Hb)
        ! accept the converged Ritz values, pushed outward by their residual
        ! (an eigenvalue lies within that distance) and by the safety margin
        nacc = 0; best_rel = HUGE(1._xp); k = 1
        DO im = 1,me
          rel = res(im)/MAX(ABS(eig(im)), TINY(1._xp))
          IF (rel .NE. rel) rel = HUGE(1._xp) ! NaN guard: never accept it
          acc(im) = rel .LE. RITZ_TOL
          IF (acc(im)) nacc = nacc + 1
          IF (rel .LT. best_rel) THEN
            best_rel = rel; k = im
          ENDIF
        ENDDO
        IF (nacc .EQ. 0) THEN
          acc(1:me) = .FALSE.; acc(k) = .TRUE.
        ENDIF
        hblk = HUGE_H
        DO im = 1,me
          hlam(im) = HUGE_H
          IF (.NOT. acc(im)) CYCLE
          rel = res(im)/MAX(ABS(eig(im)), TINY(1._xp))
          IF (rel .NE. rel) rel = 0._xp
          lam(im)  = eig(im)*(1._xp + rel)/dt_lin_safety
          hlam(im) = h_stab(lam(im))
          IF (hlam(im) .LT. hblk) THEN
            hblk = hlam(im); k = im
          ENDIF
        ENDDO
        DO im = blk_ms(ib),blk_me(ib)
          dt_lin_map(mode_iky(im),mode_ikx(im)) = hblk
        ENDDO
        IF (hblk .LT. dt_lin_loc) THEN
          dt_lin_loc = hblk
          lim_buf(1) = kyarray(mode_iky(blk_ms(ib)))
          lim_buf(2) = kxarray(mode_iky(blk_ms(ib)),mode_ikx(blk_ms(ib)))
          lim_buf(3) = REAL(eig(k),xp); lim_buf(4) = AIMAG(eig(k))
          lim_resid  = res(k)
        ENDIF
        CALL keep_extremes(me, lam(1:me), hlam(1:me), acc(1:me))
      ENDDO
      DEALLOCATE(H, m_eff)
    ENDDO
    ! restore the state exactly (the initial phi may deliberately differ from
    ! the field solve of the moments, e.g. INIT_OPT = 'phi')
    moments(:,:,:,:,:,:,1) = moments_save
    phi = phi_save; psi = psi_save; dBpar = dBpar_save
    DEALLOCATE(moments_save, phi_save, psi_save, dBpar_save)
    CALL compute_Sapj ! S(phi_0,N_0) as initialize leaves it
    DEALLOCATE(eig, lam, res, hlam, acc)
    IF (ALLOCATED(mode_iky)) DEALLOCATE(mode_iky, mode_ikx, blk_ms, blk_me, blk_bs, blk_be)
    ! global limit and the block that sets it
    CALL MPI_ALLREDUCE(dt_lin_loc, dt_lin, 1, mpi_xp_r, MPI_MIN, MPI_COMM_WORLD, ierr)
    owner = HUGE(1)
    IF (dt_lin_loc .EQ. dt_lin) owner = my_id
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, owner, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(lim_buf,   4, mpi_xp_r, owner, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST(lim_resid, 1, mpi_xp_r, owner, MPI_COMM_WORLD, ierr)
    lim_ky = lim_buf(1); lim_kx = lim_buf(2); lim_eig = CMPLX(lim_buf(3), lim_buf(4), xp)
    dt_lin_estimated = .TRUE.
    arnoldi_time = MPI_WTIME() - t0
  CONTAINS
    !> Keep the N_KEEP Ritz values of a block that matter for the omega table:
    !! the smallest stable steps, the largest |Im| and the most negative Re.
    SUBROUTINE keep_extremes(me_, lam_, hlam_, acc_)
      INTEGER, INTENT(IN) :: me_
      COMPLEX(xp), DIMENSION(me_), INTENT(IN) :: lam_
      REAL(xp),    DIMENSION(me_), INTENT(IN) :: hlam_
      LOGICAL,     DIMENSION(me_), INTENT(IN) :: acc_
      LOGICAL, DIMENSION(me_) :: taken
      INTEGER :: n_, i_, pick, crit
      REAL(xp) :: score, best
      taken = .NOT. acc_
      n_ = 0
      DO crit = 1,N_KEEP
        pick = 0; best = HUGE(1._xp)
        DO i_ = 1,me_
          IF (taken(i_)) CYCLE
          SELECT CASE (MOD(crit,3))
          CASE (1); score = hlam_(i_)            ! smallest stable step
          CASE (2); score = -ABS(AIMAG(lam_(i_))) ! largest frequency
          CASE (0); score = REAL(lam_(i_),xp)     ! strongest damping
          END SELECT
          IF (score .LT. best) THEN
            best = score; pick = i_
          ENDIF
        ENDDO
        IF (pick .EQ. 0) EXIT
        taken(pick) = .TRUE.
        n_ritz = n_ritz + 1
        ritz_list(n_ritz) = lam_(pick)
        n_ = n_ + 1
      ENDDO
    END SUBROUTINE keep_extremes
  END SUBROUTINE estimate_linear_limit

  !=============================================================================
  !                     Combined limit as a function of omega_NL
  !=============================================================================
  !> Tabulate 1/h_stab(omega) = max over the kept Ritz values lambda (and 0) of
  !! 1/h_stab(lambda +/- i omega). Linear in omega on [0, omega_tab_max], with
  !! the 1/omega asymptote of pure advection beyond. 1/h is tabulated rather
  !! than h because it is close to linear in omega.
  SUBROUTINE build_h_table
    USE mpi
    IMPLICIT NONE
    INTEGER  :: i, k, ierr
    REAL(xp) :: om, h
    COMPLEX(xp) :: shift
    IF (ALLOCATED(inv_h_tab)) DEALLOCATE(inv_h_tab)
    ALLOCATE(inv_h_tab(n_omega_tab))
    IF (dt_lin .GE. 0.5_xp*HUGE_H .OR. dt_lin .LE. 0._xp) THEN
      omega_tab_max = 1._xp
    ELSE
      omega_tab_max = 8._xp*S_imag/dt_lin
    ENDIF
    domega_tab = omega_tab_max/REAL(n_omega_tab-1,xp)
    DO i = 1,n_omega_tab
      om    = REAL(i-1,xp)*domega_tab
      shift = CMPLX(0._xp, om, xp)
      h = h_stab(shift) ! pure advection
      IF (dt_lin_estimated) THEN
        DO k = 1,n_ritz
          h = MIN(h, h_stab(ritz_list(k) + shift), h_stab(ritz_list(k) - shift))
        ENDDO
      ELSE
        ! no spectrum available: add the frequencies (conservative)
        h = MIN(h, 1._xp/(1._xp/dt_lin + om/S_imag))
      ENDIF
      IF (h .GE. 0.5_xp*HUGE_H) THEN
        inv_h_tab(i) = 0._xp
      ELSE
        inv_h_tab(i) = 1._xp/h
      ENDIF
    ENDDO
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, inv_h_tab, n_omega_tab, mpi_xp_r, MPI_MAX, MPI_COMM_WORLD, ierr)
  END SUBROUTINE build_h_table

  !> Stable step (cfl_frac = 1) for an ExB CFL frequency omega, from the table.
  FUNCTION h_of_omega(om) RESULT(h)
    IMPLICIT NONE
    REAL(xp), INTENT(IN) :: om
    REAL(xp) :: h, x, f, inv_h
    INTEGER  :: i
    IF (om .LE. 0._xp) THEN
      inv_h = inv_h_tab(1)
    ELSE IF (om .GE. omega_tab_max) THEN
      inv_h = inv_h_tab(n_omega_tab)*om/omega_tab_max
    ELSE
      x = om/domega_tab
      i = MIN(INT(x) + 1, n_omega_tab - 1)
      f = x - REAL(i-1,xp)
      inv_h = (1._xp - f)*inv_h_tab(i) + f*inv_h_tab(i+1)
    ENDIF
    IF (inv_h .LE. 0._xp) THEN
      h = HUGE_H
    ELSE
      h = 1._xp/inv_h
    ENDIF
  END FUNCTION h_of_omega

  !=============================================================================
  !                         ExB CFL frequency of a state
  !=============================================================================
  !> omega_nl = max over z and the (x,y) plane of kx_eff |d_y F| + ky_eff |d_x F|
  !! for the potentials the Poisson bracket advects with: phi (plus the
  !! background zonal mode), tau/q dBpar, and psi scaled by the parallel
  !! velocity bound of the flutter bracket. The gyroaverage kernels are bounded
  !! by one, so the bare potentials give an upper bound. Zero in linear runs.
  SUBROUTINE compute_omega_nl
    USE fields,  ONLY: phi, psi, dBpar
    USE grid,    ONLY: local_nky, local_nkx, local_nz, ngz_o2, total_nkx, Ny, &
                       kxarray, kyarray, AA_x, AA_y, iky0, contains_kx0, contains_ky0, deltakx
    USE model,   ONLY: LINEARITY, EM, BPAR, ZFrate, ikxZF, ZF_ONLY
    USE species, ONLY: tau_q
    USE fourier, ONLY: grad_xy_max
    USE mpi
    IMPLICIT NONE
    COMPLEX(xp), DIMENSION(local_nky,local_nkx) :: F
    REAL(xp), DIMENSION(2) :: keff
    REAL(xp) :: kx_eff, ky_eff, w, wz, wloc, tauq_max, zf_amp
    INTEGER  :: iky, ikx, iz, izi, ierr, ikxExBp, ikxExBn
    omega_nl = 0._xp
    IF (TRIM(LINEARITY) .NE. 'nonlinear') RETURN
    ! largest wavenumbers the bracket keeps (kx moves under ExB shear)
    keff = 0._xp
    DO ikx = 1,local_nkx
      DO iky = 1,local_nky
        IF (AA_x(ikx)*AA_y(iky) .LE. 0._xp) CYCLE
        keff(1) = MAX(keff(1), ABS(kxarray(iky,ikx)))
        keff(2) = MAX(keff(2), ABS(kyarray(iky)))
      ENDDO
    ENDDO
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, keff, 2, mpi_xp_r, MPI_MAX, MPI_COMM_WORLD, ierr)
    kx_eff = keff(1); ky_eff = keff(2)
    tauq_max = 0._xp
    IF (BPAR) tauq_max = MAXVAL(ABS(tau_q))
    wloc = 0._xp
    DO iz = 1,local_nz
      izi = iz + ngz_o2
      F = phi(:,:,izi)
      IF (ABS(ZFrate) .GT. 0._xp) THEN ! background zonal mode, as in compute_nonlinear
        ikxExBp = ikxZF + 1
        ikxExBn = total_nkx - (ikxExBp - 2)
        zf_amp  = 8._xp*ZFrate*REAL(ikxZF,xp)/deltakx
        IF (ZF_ONLY) F = (0._xp, 0._xp)
        IF (contains_kx0 .AND. contains_ky0) THEN
          F(iky0,ikxExBp) = F(iky0,ikxExBp) + zf_amp
          F(iky0,ikxExBn) = F(iky0,ikxExBn) + zf_amp
        ENDIF
      ENDIF
      CALL grad_xy_max(kyarray, kxarray, AA_y, AA_x, local_nky, total_nkx, Ny, F, kx_eff, ky_eff, w)
      wz = w
      IF (BPAR) THEN
        F = dBpar(:,:,izi)
        CALL grad_xy_max(kyarray, kxarray, AA_y, AA_x, local_nky, total_nkx, Ny, F, kx_eff, ky_eff, w)
        wz = wz + tauq_max*w
      ENDIF
      IF (EM) THEN
        F = psi(:,:,izi)
        CALL grad_xy_max(kyarray, kxarray, AA_y, AA_x, local_nky, total_nkx, Ny, F, kx_eff, ky_eff, w)
        wz = wz + flutter_factor*w
      ENDIF
      wloc = MAX(wloc, wz)
    ENDDO
    CALL MPI_ALLREDUCE(wloc, omega_nl, 1, mpi_xp_r, MPI_MAX, MPI_COMM_WORLD, ierr)
  END SUBROUTINE compute_omega_nl

END MODULE adaptive_dt
