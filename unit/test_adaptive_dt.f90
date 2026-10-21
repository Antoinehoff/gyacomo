PROGRAM test_adaptive_dt
  ! Verifies the building blocks of adaptive_dt_mod and, on the real RHS of
  ! the code, the two measurements it relies on.
  !
  ! The mode is taken from the command line because readinputs and auxval can
  ! only run once per process; "make check" runs this binary once per mode.
  !
  !   usage: test_adaptive_dt [hyperdiff|cfl]
  !
  !   hyperdiff: a linear z-pinch deck where every (kx,ky) block of the linear
  !              operator is 1x1 with the eigenvalue of the x/y hyperdiffusion,
  !              -mu_x (kx/kx_max)^N_HD - mu_y (ky/ky_max)^N_HD, so that the
  !              Arnoldi estimate must be exact and the state must be restored.
  !   cfl:       a nonlinear deck with a single-mode phi of known gradient, so
  !              that omega_nl pins the c2r normalization and the padding rows.
  !
  ! Both modes also run the scheme-only checks (stability polynomial, extents,
  ! Hessenberg eigenvalues, Arnoldi kernel on an injected operator, controller).
  USE prec_const,       ONLY: xp, imagu, PI, SQRT3
  USE basic,            ONLY: lu_in, change_dt, dt
  USE parallel,         ONLY: ppinit, comm_pz, my_id
  USE grid,             ONLY: local_nky, local_nkx, local_nz, total_nky, total_nkx, ngz_o2, &
                              kxarray, kyarray, kx_max, ky_max, diff_kx_coeff, diff_ky_coeff, &
                              deltakx, deltaky, AA_x, AA_y, ip0, ij0, ngp_o2, ngj_o2, Nx, Ny, iky0, ikx0
  USE model,            ONLY: mu_x, mu_y, N_HD
  USE fields,           ONLY: moments, phi, psi, dBpar
  USE time_integration, ONLY: set_updatetlevel, numerical_scheme, ntimelevel
  USE adaptive_dt
  USE test_harness
  IMPLICIT NONE

  CHARACTER(len=32) :: mode
  INTEGER  :: nargs, iky, ikx, k, n, info, i
  REAL(xp) :: tol, want, got, h, x, y, dt_new, A_amp, dx, dy
  COMPLEX(xp) :: z, R, lam
  COMPLEX(xp), DIMENSION(5,5) :: Cm
  COMPLEX(xp), DIMENSION(5)   :: ev
  COMPLEX(xp), DIMENSION(:),   ALLOCATABLE :: eig
  REAL(xp),    DIMENSION(:),   ALLOCATABLE :: res
  COMPLEX(xp), DIMENSION(:,:,:), ALLOCATABLE :: Hm
  COMPLEX(xp), DIMENSION(:,:),   ALLOCATABLE :: Hb
  INTEGER,     DIMENSION(:),   ALLOCATABLE :: m_eff
  COMPLEX(xp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: mom_before
  COMPLEX(xp), DIMENSION(:,:,:),       ALLOCATABLE :: phi_before
  ! injected diagonal operator for the Arnoldi kernel test
  INTEGER, PARAMETER :: NDIAG = 40
  COMPLEX(xp), DIMENSION(NDIAG) :: diag_op
  LOGICAL :: ok

  CALL ppinit
  nargs = COMMAND_ARGUMENT_COUNT()
  IF (nargs .GE. 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, mode)
  ELSE
    mode = 'hyperdiff'
  ENDIF
  SELECT CASE (TRIM(mode))
  CASE ('hyperdiff', 'cfl')
  CASE DEFAULT
    WRITE(*,'(A)') 'unknown mode "'//TRIM(mode)//'"'
    CALL start_suite('adaptive_dt ['//TRIM(mode)//']')
    CALL check('mode is known to the test', .FALSE.)
    CALL finish_suite
    STOP
  END SELECT

  CALL start_suite('adaptive_dt ['//TRIM(mode)//']')
  tol = tol_of(1.0e4_xp)

  CALL write_input_deck
  CALL readinputs
  CALL auxval
  CALL set_updatetlevel(1)
  CALL change_dt(0.01_xp)

  !------------------------------------------------- scheme characterization
  CALL section('stability polynomial of the tableau (RK4)')
  CALL check('deck selected RK4', TRIM(numerical_scheme) .EQ. 'RK4')
  ok = .TRUE.
  DO k = 1,5
    z = CMPLX(-0.7_xp*k, 0.4_xp*k, xp)
    R = stability_function(z)
    want = ABS(1._xp + z + z**2/2._xp + z**3/6._xp + z**4/24._xp)
    IF (ABS(ABS(R) - want) .GT. tol*MAX(1._xp, want)) ok = .FALSE.
  ENDDO
  CALL check('R(z) is the degree 4 Taylor polynomial', ok)

  CALL setup_stability
  CALL check_real('S_real of RK4', S_real, 2.7852935634_xp, 1.0e-6_xp)
  CALL check_real('S_imag of RK4', S_imag, 2._xp*SQRT(2._xp), 1.0e-6_xp)

  CALL section('largest stable step per eigenvalue')
  CALL check_real('h_stab(-1) = S_real', h_stab((-1._xp, 0._xp)), S_real, 1.0e-6_xp)
  CALL check_real('h_stab(i)  = S_imag', h_stab((0._xp, 1._xp)),  S_imag, 1.0e-6_xp)
  CALL check_real('h_stab(-2) = S_real/2', h_stab((-2._xp, 0._xp)), 0.5_xp*S_real, 1.0e-6_xp)
  CALL check_real('unstable eigenvalue projected on the imaginary axis', &
                  h_stab((0.5_xp, 1._xp)), S_imag, 1.0e-6_xp)
  CALL check('h_stab(0) means no limit', h_stab((0._xp, 0._xp)) .GT. 1.0e20_xp)

  !------------------------------------------------- small linear algebra
  CALL section('Hessenberg eigenvalues')
  ! companion matrix of (z-1)(z-2)(z-3)(z-(1+2i))(z-(1-2i))
  !   = z^5 - 8 z^4 + 28 z^3 - 58 z^2 + 67 z - 30
  Cm = (0._xp, 0._xp)
  Cm(1,1) = (8._xp,0._xp); Cm(1,2) = (-28._xp,0._xp); Cm(1,3) = (58._xp,0._xp)
  Cm(1,4) = (-67._xp,0._xp); Cm(1,5) = (30._xp,0._xp)
  DO k = 2,5
    Cm(k,k-1) = (1._xp, 0._xp)
  ENDDO
  CALL hessenberg_eig(5, Cm, ev, info)
  CALL check_int('QR converged on the companion matrix', info, 0)
  ok = found(ev, (1._xp,0._xp)) .AND. found(ev, (2._xp,0._xp)) .AND. found(ev, (3._xp,0._xp)) &
       .AND. found(ev, (1._xp,2._xp)) .AND. found(ev, (1._xp,-2._xp))
  CALL check('roots 1, 2, 3, 1+2i, 1-2i recovered', ok)
  ! upper triangular: eigenvalues are the diagonal
  Cm = (0._xp, 0._xp)
  DO k = 1,5
    Cm(k,k) = CMPLX(REAL(k,xp), -0.3_xp*k, xp)
    DO i = k+1,5
      Cm(k,i) = CMPLX(0.1_xp*i, 0.2_xp*k, xp)
    ENDDO
  ENDDO
  CALL hessenberg_eig(5, Cm, ev, info)
  ok = .TRUE.
  DO k = 1,5
    ok = ok .AND. found(ev, CMPLX(REAL(k,xp), -0.3_xp*k, xp))
  ENDDO
  CALL check('triangular matrix: diagonal recovered', ok)
  CALL check_real('Hermite ladder radius P=2 is sqrt(3)', hermite_ladder_radius(2), SQRT3, 1.0e-10_xp)
  CALL check_real('Hermite ladder radius P=1 is 1',       hermite_ladder_radius(1), 1._xp, 1.0e-10_xp)

  !------------------------------------------------- Arnoldi kernel
  CALL section('Arnoldi on an injected diagonal operator')
  DO k = 1,NDIAG
    ! eigenvalues on an arc in the left half plane, well separated extremes
    x = -0.05_xp*k
    y = 3._xp*SIN(PI*REAL(k,xp)/REAL(NDIAG+1,xp))
    diag_op(k) = CMPLX(x, y, xp)
  ENDDO
  diag_op(NDIAG) = (-6._xp, 0._xp) ! most damped, isolated (gap 4 from the arc)
  diag_op(1)     = (-0.1_xp, 5._xp) ! fastest, isolated (gap 2 from the arc)
  ! full dimension: exact spectrum
  ALLOCATE(Hm(NDIAG+1,NDIAG,1), m_eff(1), eig(NDIAG), res(NDIAG))
  CALL arnoldi_blocks(diag_matvec, NDIAG, 1, (/1/), (/NDIAG/), NDIAG, comm_pz, 7, Hm, m_eff)
  CALL check_int('full Krylov space reached (m_eff = n)', m_eff(1), NDIAG)
  ALLOCATE(Hb(NDIAG,NDIAG)); Hb = Hm(1:NDIAG,1:NDIAG,1)
  CALL hessenberg_eig(NDIAG, Hb, eig, info)
  ok = .TRUE.
  DO k = 1,NDIAG
    ok = ok .AND. found(eig, diag_op(k), 1.0e-8_xp)
  ENDDO
  CALL check('m = n recovers every eigenvalue', ok)
  CALL ritz_residuals(NDIAG, Hb, ABS(Hm(NDIAG+1,NDIAG,1)), eig, res)
  CALL check('residuals vanish on an invariant space', MAXVAL(res) .LT. 1.0e-8_xp)
  DEALLOCATE(Hm, m_eff, eig, res, Hb)
  ! small Krylov space: the isolated extremes converge first
  n = 16
  ALLOCATE(Hm(n+1,n,1), m_eff(1), eig(n), res(n))
  CALL arnoldi_blocks(diag_matvec, NDIAG, 1, (/1/), (/NDIAG/), n, comm_pz, 7, Hm, m_eff)
  ALLOCATE(Hb(n,n)); Hb = Hm(1:n,1:n,1)
  CALL hessenberg_eig(n, Hb, eig, info)
  CALL ritz_residuals(n, Hb, ABS(Hm(n+1,n,1)), eig, res)
  CALL check('m < n: most damped eigenvalue found to 1e-5',  found(eig, (-6._xp, 0._xp),  1.0e-5_xp))
  CALL check('m < n: fastest eigenvalue found to 1e-5',      found(eig, (-0.1_xp, 5._xp), 1.0e-5_xp))
  ok = .TRUE.
  DO k = 1,n
    IF (ABS(eig(k) - (-6._xp,0._xp)) .LT. 1.0e-4_xp .OR. ABS(eig(k) - (-0.1_xp,5._xp)) .LT. 1.0e-4_xp) THEN
      IF (res(k) .GT. 1.0e-3_xp) ok = .FALSE.
    ENDIF
  ENDDO
  CALL check('their residuals are small', ok)
  DEALLOCATE(Hm, m_eff, eig, res, Hb)

  !------------------------------------------------- controller
  CALL section('step size controller')
  ! defaults: dt_growth = 1.2, dt_hyst = 0.1, dt_min = 1e-6, dt_max = 1e6
  CALL check_real('decrease is immediate',     next_dt(1._xp, 0.3_xp, 0._xp, 100._xp), 0.3_xp)
  CALL check_real('growth capped at dt_growth', next_dt(1._xp, 5._xp, 0._xp, 100._xp), 1.2_xp)
  CALL check_real('dead zone below dt_hyst',   next_dt(1._xp, 1.05_xp, 0._xp, 100._xp), 1._xp)
  CALL check_real('small increase taken as is', next_dt(1._xp, 1.15_xp, 0._xp, 100._xp), 1.15_xp)
  CALL check_real('clipped to dt_max',         next_dt(1.0e6_xp, 5.0e6_xp, 0._xp, 1.0e9_xp), 1.0e6_xp)
  dt_new = next_dt(1._xp, 1._xp, 99.5_xp, 100._xp)
  CALL check('lands on tmax', (dt_new .GE. 0.5_xp) .AND. (99.5_xp + dt_new .GE. 100._xp) &
                              .AND. (dt_new .LT. 0.5_xp*(1._xp + 1.0e-6_xp)))

  !------------------------------------------------- on the real RHS
  SELECT CASE (TRIM(mode))
  CASE ('hyperdiff')
    CALL section('linear limit on a pure hyperdiffusion deck')
    CALL check_int('single rank: local_nky is the total', local_nky, total_nky)
    CALL check_int('single rank: local_nkx is the total', local_nkx, total_nkx)
    ! something to restore: a nonzero state with phi deliberately inconsistent
    moments = (0._xp, 0._xp)
    moments(1,ip0,ij0,:,:,1+ngz_o2,1) = (0.3_xp, -0.1_xp)
    phi = (0.7_xp, 0.2_xp); psi = (0._xp, 0._xp); dBpar = (0._xp, 0._xp)
    ALLOCATE(mom_before(SIZE(moments,1),SIZE(moments,2),SIZE(moments,3),SIZE(moments,4),SIZE(moments,5),SIZE(moments,6)))
    ALLOCATE(phi_before(SIZE(phi,1),SIZE(phi,2),SIZE(phi,3)))
    mom_before = moments(:,:,:,:,:,:,1); phi_before = phi
    CALL estimate_linear_limit
    CALL check('blocks are (kx,ky) modes in a z-pinch', .NOT. row_blocks)
    CALL check('state restored: moments', ALL(moments(:,:,:,:,:,:,1) .EQ. mom_before))
    CALL check('state restored: phi',     ALL(phi .EQ. phi_before))
    ! every block is 1x1: the Arnoldi breaks down at the first step and the
    ! Ritz value is the exact eigenvalue, so the per mode limit is
    ! dt_lin_safety * h_stab(lambda) = dt_lin_safety * S_real/|lambda|
    ok = .TRUE.
    DO ikx = 1,local_nkx
      DO iky = 1,local_nky
        lam  = CMPLX(-mu_x*diff_kx_coeff*kxarray(iky,ikx)**N_HD - mu_y*diff_ky_coeff*kyarray(iky)**N_HD, 0._xp, xp)
        IF (ABS(lam) .LT. 1.0e-12_xp) THEN
          want = 1.0e30_xp
        ELSE
          want = dt_lin_safety*S_real/ABS(lam)
        ENDIF
        got = dt_lin_map(iky,ikx)
        IF (want .GT. 1.0e20_xp) THEN
          IF (got .LT. 1.0e20_xp) ok = .FALSE.
        ELSE
          IF (ABS(got - want) .GT. 1.0e-6_xp*want) ok = .FALSE.
        ENDIF
      ENDDO
    ENDDO
    CALL check('per mode limit equals S_real/|lambda| (with the safety margin)', ok)
    want = dt_lin_safety*S_real/(mu_x + mu_y) ! corner mode kx = kx_max, ky = ky_max
    CALL check_real('global dt_lin set by the corner mode', dt_lin, want, 1.0e-6_xp)
    CALL check_real('limiting eigenvalue real part', REAL(lim_eig,xp), -(mu_x + mu_y), 1.0e-6_xp)
    CALL check('limiting residual is zero (exact 1x1 blocks)', lim_resid .LT. 1.0e-10_xp)
    ! the omega table: no advection reproduces dt_lin, large omega the pure CFL
    CALL build_h_table
    CALL check_real('h(omega=0) = dt_lin', h_of_omega(0._xp), dt_lin, 1.0e-6_xp)
    x = 50._xp*S_imag/dt_lin
    got = h_of_omega(x)
    CALL check('h(large omega) tends to S_imag/omega', ABS(got - S_imag/x) .LT. 0.1_xp*S_imag/x)
    CALL check('h(omega) is monotone', h_of_omega(0.5_xp*x) .GT. got .AND. h_of_omega(0.5_xp*x) .LT. dt_lin)
  CASE ('cfl')
    CALL section('ExB CFL frequency of a single mode potential')
    ! phi = A cos(k0 x) with k0 the first kx: the (ky=0, kx=+-k0) pair with
    ! amplitude A/2 each. d_x phi has amplitude A k0 and d_y phi vanishes, so
    ! omega_nl = ky_eff * A * k0 with ky_eff the largest anti-aliased ky.
    A_amp = 0.8_xp
    phi = (0._xp, 0._xp); psi = (0._xp, 0._xp); dBpar = (0._xp, 0._xp)
    phi(iky0,ikx0+1,:)       = CMPLX(0.5_xp*A_amp, 0._xp, xp)
    phi(iky0,total_nkx,:)    = CMPLX(0.5_xp*A_amp, 0._xp, xp)
    CALL compute_omega_nl
    want = 0._xp
    DO iky = 1,local_nky
      IF (AA_y(iky) .GT. 0._xp) want = MAX(want, ABS(kyarray(iky)))
    ENDDO
    want = want*A_amp*deltakx
    CALL check_real('omega_nl = ky_eff A k0', omega_nl, want, 1.0e-10_xp)
    ! a mode beyond the 2/3 cutoff must not count
    phi = (0._xp, 0._xp)
    phi(iky0,Nx/2+1,:) = CMPLX(A_amp, 0._xp, xp) ! kx = kx_max, filtered
    CALL compute_omega_nl
    CALL check_real('anti-aliased modes do not advect', omega_nl, 0._xp, 1.0e-12_xp)
    ! y dependence: phi = A cos(k0y y): d_y phi amplitude A k0y, omega = kx_eff A k0y
    phi = (0._xp, 0._xp)
    phi(iky0+1,ikx0,:) = CMPLX(0.5_xp*A_amp, 0._xp, xp) ! ky = dky, its conjugate is implicit
    CALL compute_omega_nl
    want = 0._xp
    DO ikx = 1,local_nkx
      DO iky = 1,local_nky
        IF (AA_x(ikx)*AA_y(iky) .GT. 0._xp) want = MAX(want, ABS(kxarray(iky,ikx)))
      ENDDO
    ENDDO
    want = want*A_amp*deltaky
    CALL check_real('omega_nl = kx_eff A k0y', omega_nl, want, 1.0e-10_xp)
  END SELECT

  CALL finish_suite

CONTAINS

  LOGICAL FUNCTION found(list, val, tolerance)
    COMPLEX(xp), DIMENSION(:), INTENT(IN) :: list
    COMPLEX(xp), INTENT(IN) :: val
    REAL(xp), OPTIONAL, INTENT(IN) :: tolerance
    REAL(xp) :: t_
    t_ = 1.0e-8_xp
    IF (PRESENT(tolerance)) t_ = tolerance
    found = ANY(ABS(list - val) .LT. t_*MAX(1._xp, ABS(val)))
  END FUNCTION found

  SUBROUTINE diag_matvec(v, w)
    COMPLEX(xp), DIMENSION(:), INTENT(IN)  :: v
    COMPLEX(xp), DIMENSION(:), INTENT(OUT) :: w
    w = diag_op*v
  END SUBROUTINE diag_matvec

  !> Feed the namelists through a scratch file so nothing is left on disk.
  SUBROUTINE write_input_deck
    OPEN(UNIT=lu_in, STATUS='SCRATCH', ACTION='READWRITE', FORM='FORMATTED')
    WRITE(lu_in,'(A)')       '&GRID'
    WRITE(lu_in,'(A)')       '  pmax = 0'
    WRITE(lu_in,'(A)')       '  jmax = 0'
    WRITE(lu_in,'(A)')       '  Nx   = 16'
    WRITE(lu_in,'(A)')       '  Lx   = 40.0'
    WRITE(lu_in,'(A)')       '  Ny   = 12'
    WRITE(lu_in,'(A)')       '  Ly   = 30.0'
    WRITE(lu_in,'(A)')       '  Nz   = 1'
    WRITE(lu_in,'(A)')       '  SG   = .false.'
    WRITE(lu_in,'(A)')       '  Nexc = 1'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&GEOMETRY'
    WRITE(lu_in,'(A)')       "  geom = 'z-pinch'"
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&DIAGNOSTICS'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&MODEL'
    IF (TRIM(mode) .EQ. 'cfl') THEN
      WRITE(lu_in,'(A)')     "  LINEARITY = 'nonlinear'"
    ELSE
      WRITE(lu_in,'(A)')     "  LINEARITY = 'linear'"
    ENDIF
    WRITE(lu_in,'(A)')       '  Na      = 1'
    WRITE(lu_in,'(A)')       '  ADIAB_E = .true.'
    WRITE(lu_in,'(A)')       '  mu_x    = 0.7'
    WRITE(lu_in,'(A)')       '  mu_y    = 1.3'
    WRITE(lu_in,'(A)')       '  N_HD    = 4'
    WRITE(lu_in,'(A)')       '  nu      = 0.0'
    WRITE(lu_in,'(A)')       '  k_cB    = 0.0'
    WRITE(lu_in,'(A)')       '  k_gB    = 0.0'
    WRITE(lu_in,'(A)')       '  beta    = 0.0'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&CLOSURE'
    WRITE(lu_in,'(A)')       "  hierarchy_closure = 'truncation'"
    WRITE(lu_in,'(A)')       "  nonlinear_closure = 'truncation'"
    WRITE(lu_in,'(A)')       '  nmax = 0'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&SPECIES'
    WRITE(lu_in,'(A)')       "  name_  = 'ions'"
    WRITE(lu_in,'(A)')       '  tau_   = 1.0'
    WRITE(lu_in,'(A)')       '  sigma_ = 1.0'
    WRITE(lu_in,'(A)')       '  q_     = 1.0'
    WRITE(lu_in,'(A)')       '  k_N_   = 0.0'
    WRITE(lu_in,'(A)')       '  k_T_   = 0.0'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&COLLISION'
    WRITE(lu_in,'(A)')       "  collision_model = 'DG'"
    WRITE(lu_in,'(A)')       '  GK_CO           = .false.'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&INITIAL'
    WRITE(lu_in,'(A)')       "  INIT_OPT = 'phi'"
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&TIME_INTEGRATION'
    WRITE(lu_in,'(A)')       "  numerical_scheme = 'RK4'"
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&ADAPTIVE_DT'
    WRITE(lu_in,'(A)')       '  use_adaptive_dt = .true.'
    WRITE(lu_in,'(A)')       '  krylov_dim      = 8'
    WRITE(lu_in,'(A)')       '/'
    REWIND(lu_in)
  END SUBROUTINE write_input_deck

END PROGRAM test_adaptive_dt
