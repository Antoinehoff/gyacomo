PROGRAM test_bpar
  ! Exercises the parallel magnetic fluctuation (delta B_par) machinery:
  ! the K'_j = K_j - K_{j-1} kernels, the coupled quasineutrality / pressure
  ! balance operator, the field solve, the non adiabatic moments and the
  ! gradient drive of the moment hierarchy.
  !
  ! Everything is measured against closed forms written in this file (the
  ! kernels y^j e^-y/j!, Bessel functions from the intrinsics, the 2x2 operator
  ! built by hand) rather than against the code's own arrays, so that a sign
  ! or a factor slipped into numerics_mod cannot pass by being self-consistent.
  ! The one exception is the field equation residual, which is evaluated with
  ! the code's kernels on purpose: it pins the solve to the non adiabatic
  ! moments it has to be consistent with, whatever the kernel model.
  !
  ! The mode is taken from the command line because readinputs can only run
  ! once per process (allocate_array does not deallocate), so "make check" runs
  ! this binary once per mode:
  !
  !   on     : beta > 0, BPAR = .true.
  !   off    : beta > 0, BPAR = .false., the A_par-only path must be untouched
  !   nobeta : beta = 0, BPAR = .true. is demoted to off at input time
  !
  ! Single rank, z-pinch with a single z plane: every local size is the total,
  ! hatB = 1, and the parallel terms of the hierarchy vanish identically.
  USE prec_const, ONLY: xp, imagu, SQRT2
  USE basic,      ONLY: lu_in
  USE parallel,   ONLY: ppinit
  USE grid,       ONLY: total_nkx, total_nky, local_nkx, local_nky, local_na, local_np, local_nj, &
                        ngj, ngj_o2, ngp_o2, ngz_o2, jarray, parray, kxarray, kyarray, kp2array, &
                        ieven, ip0, ij0, ikx0, iky0, jmax
  USE fields,     ONLY: moments, phi, psi, dBpar
  USE array,      ONLY: kernel, kernel_dB, inv_poisson_op, bp_phi_rho, bp_xcoupl, bp_dB_pp, &
                        moments_rhs
  USE model,      ONLY: BPAR, beta
  USE species,    ONLY: q, tau, k_N, k_T, tau_q, sigma2_tau_o2
  USE geometry,   ONLY: hatB
  USE time_integration, ONLY: updatetlevel, set_updatetlevel
  USE processing, ONLY: compute_nadiab_moments, compute_gradients_z, compute_interp_z
  USE nonlinear,  ONLY: compute_Sapj
  USE moments_eq_rhs, ONLY: compute_moments_eq_rhs
  USE test_harness
  IMPLICIT NONE

  ! jmax = 16 keeps the truncated kernel sums within round off of the Bessel
  ! closed forms for the y <= 2 modes used below (tail ~ y^j/j!).
  INTEGER,  PARAMETER :: T_PMAX = 4, T_JMAX = 16
  INTEGER,  PARAMETER :: T_NX = 8, T_NY = 8, T_NZ = 1
  REAL(xp), PARAMETER :: T_LX = 40._xp, T_LY = 40._xp
  REAL(xp), PARAMETER :: T_BETA = 0.05_xp
  COMPLEX(xp), PARAMETER :: B_AMP = (0.3_xp, -0.2_xp) ! dBpar amplitude for the drive check

  CHARACTER(len=32) :: mode
  REAL(xp) :: tol, beta_in, y, x, lhs, rhs, Kc(-1:T_JMAX+1), Lag(0:T_JMAX+1), lam
  REAL(xp) :: P_, C_, M_, D_, sKK, sKKp, sKpKp, worst, scale
  COMPLEX(xp) :: rho, pperp, hj, R_qn, R_B, got, want
  INTEGER  :: nargs, ia, ij, iji, ip, ipi, iky, ikx, izi, p_int, j_int, nmode
  LOGICAL  :: bpar_wanted, ok

  CALL ppinit
  nargs = COMMAND_ARGUMENT_COUNT()
  IF (nargs .GE. 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, mode)
  ELSE
    mode = 'on'
  ENDIF
  SELECT CASE (TRIM(mode))
  CASE ('on');     bpar_wanted = .TRUE.;  beta_in = T_BETA
  CASE ('off');    bpar_wanted = .FALSE.; beta_in = T_BETA
  CASE ('nobeta'); bpar_wanted = .TRUE.;  beta_in = 0._xp
  CASE DEFAULT
    WRITE(*,'(A)') 'unknown mode "'//TRIM(mode)//'"'
    CALL start_suite('bpar ['//TRIM(mode)//']')
    CALL check('mode is known to the test', .FALSE.)
    CALL finish_suite
    STOP
  END SELECT

  CALL start_suite('bpar ['//TRIM(mode)//']')
  tol = tol_of(1.0e4_xp)
  izi = 1 + ngz_o2

  CALL write_input_deck
  CALL readinputs
  CALL auxval
  CALL set_updatetlevel(1) ! the time level the solver reads, set by initial_mod in a run

  !--------------------------------------------------------------- Setup
  CALL section('the run was configured the way the test assumes')
  CALL check_int('single rank: local_nky is the total', local_nky, total_nky)
  CALL check_int('single rank: local_nkx is the total', local_nkx, total_nkx)
  CALL check_int('two kinetic species', local_na, 2)
  CALL check_int('jmax as requested', jmax, T_JMAX)
  CALL check_real('beta as requested', beta, beta_in, tol)
  CALL check_real('z-pinch has hatB = 1', hatB(izi,ieven), 1._xp, tol)
  IF (TRIM(mode) .EQ. 'nobeta') THEN
    CALL check('BPAR without beta is demoted to off', .NOT. BPAR)
  ELSE
    CALL check('BPAR read back', BPAR .EQV. bpar_wanted)
  ENDIF

  !------------------------------------------------------ kernel_dB
  ! Filled for every mode (it is cheap and keeps the code path single), so it
  ! can be pinned down regardless of the mode.
  CALL section('kernel_dB is K_j - K_{j-1} of the closed form kernel')
  ok = .TRUE.; worst = 0._xp
  DO ia = 1,local_na
    DO ikx = 1,local_nkx
      DO iky = 1,local_nky
        y = sigma2_tau_o2(ia)*kp2array(iky,ikx,izi,ieven)
        CALL closed_form_kernels(y, Kc)
        DO ij = 1,local_nj+ngj
          j_int = jarray(ij)
          IF (j_int .LT. 0) THEN
            want = 0._xp
          ELSE
            want = Kc(j_int) - Kc(j_int-1)
          ENDIF
          worst = MAX(worst, ABS(kernel_dB(ia,ij,iky,ikx,izi,ieven) - REAL(want,xp)))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  CALL check('max deviation from K_j - K_{j-1} is round off', worst .LT. tol)
  CALL check('the lower ghost (j = -1) is exactly zero', &
             ALL(kernel_dB(:,1,:,:,izi,ieven) .EQ. 0._xp))
  ! The upper ghost is reached by the j = jmax drive (xphijp1 K'_{jmax+1}).
  ia = 1; iky = 2; ikx = 2
  y = sigma2_tau_o2(ia)*kp2array(iky,ikx,izi,ieven)
  CALL closed_form_kernels(y, Kc)
  CALL check_real('the upper ghost is K_{jmax+1} - K_{jmax}', &
                  kernel_dB(ia,local_nj+ngj,iky,ikx,izi,ieven), Kc(T_JMAX+1)-Kc(T_JMAX), tol)

  CALL section('sum_j kernel_dB(j) L_j(x) is x 2J1(l)/l, l = 2 sqrt(x y)')
  ! Independent of the kernel algebra: the Laguerre series of the gyroaveraged
  ! mu delta B_par term is compared with the Bessel function it stands for.
  ! The same check on kernel pins the convention J0(l) = sum_j K_j L_j the
  ! identity rests on.
  nmode = 0
  DO ia = 1,local_na
    DO ikx = 1,local_nkx
      DO iky = 1,local_nky
        y = sigma2_tau_o2(ia)*kp2array(iky,ikx,izi,ieven)
        IF ((y .GT. 2._xp) .OR. (y .LT. 1.0e-3_xp)) CYCLE
        nmode = nmode + 1
        DO ip = 1,2
          x = MERGE(0.7_xp, 1.9_xp, ip .EQ. 1)
          CALL laguerre(x, Lag)
          lam = 2._xp*SQRT(x*y)
          lhs = 0._xp; rhs = 0._xp
          DO ij = 1,local_nj
            iji = ij + ngj_o2
            lhs = lhs + kernel_dB(ia,iji,iky,ikx,izi,ieven)*Lag(jarray(iji))
            rhs = rhs + kernel   (ia,iji,iky,ikx,izi,ieven)*Lag(jarray(iji))
          ENDDO
          CALL check_real('J0 series', rhs, BESSEL_J0(lam), tol)
          CALL check_real('x 2J1/l series', lhs, x*2._xp*BESSEL_J1(lam)/lam, tol)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  CALL check('some modes were in the convergence window', nmode .GT. 4)

  !------------------------------------------------------ operator
  IF (BPAR) THEN
    CALL section('the coupled operator matches the hand built 2x2 inverse')
    worst = 0._xp
    DO ikx = 1,local_nkx
      DO iky = 1,local_nky
        IF ((ikx .EQ. ikx0) .AND. (iky .EQ. iky0)) CYCLE
        P_ = 0._xp; C_ = 0._xp; M_ = 1._xp ! hatB = 1
        DO ia = 1,local_na
          y = sigma2_tau_o2(ia)*kp2array(iky,ikx,izi,ieven)
          CALL closed_form_kernels(y, Kc)
          sKK = 0._xp; sKKp = 0._xp; sKpKp = 0._xp
          DO j_int = 0,T_JMAX
            sKK   = sKK   + Kc(j_int)**2
            sKKp  = sKKp  + Kc(j_int)*(Kc(j_int)-Kc(j_int-1))
            sKpKp = sKpKp + (Kc(j_int)-Kc(j_int-1))**2
          ENDDO
          P_ = P_ + q(ia)**2/tau(ia)*(1._xp - sKK)
          C_ = C_ + q(ia)*sKKp
          M_ = M_ + 0.5_xp*beta*tau(ia)*sKpKp
        ENDDO
        D_ = P_*M_ + 0.5_xp*beta*C_**2
        worst = MAX(worst, ABS(bp_phi_rho(iky,ikx,1) -  M_/D_), &
                           ABS(bp_xcoupl (iky,ikx,1) + 0.5_xp*beta*C_/D_), &
                           ABS(bp_dB_pp  (iky,ikx,1) + 0.5_xp*beta*P_/D_))
        IF (D_ .LE. 0._xp) worst = HUGE(1._xp)
      ENDDO
    ENDDO
    CALL check('max deviation of the three coefficients is round off', worst .LT. tol)
    CALL check('the origin is singular and zeroed', &
               (bp_phi_rho(iky0,ikx0,1) .EQ. 0._xp) .AND. (bp_xcoupl(iky0,ikx0,1) .EQ. 0._xp) &
               .AND. (bp_dB_pp(iky0,ikx0,1) .EQ. 0._xp))
    ! Gamma_0 - Gamma_1 closed forms of the couplings, at large jmax
    CALL section('the couplings are the Gamma_0 - Gamma_1 of the literature')
    ia = 1
    DO ikx = 2,3
      DO iky = 2,3
        y = sigma2_tau_o2(ia)*kp2array(iky,ikx,izi,ieven)
        IF (y .GT. 2._xp) CYCLE
        CALL closed_form_kernels(y, Kc)
        sKKp = 0._xp; sKpKp = 0._xp
        DO j_int = 0,T_JMAX
          sKKp  = sKKp  + Kc(j_int)*(Kc(j_int)-Kc(j_int-1))
          sKpKp = sKpKp + (Kc(j_int)-Kc(j_int-1))**2
        ENDDO
        CALL check_real('sum K K''  = e^-2y (I0 - I1)(2y)', sKKp, &
                        EXP(-2._xp*y)*(bessel_i(0,2._xp*y)-bessel_i(1,2._xp*y)), tol)
        CALL check_real('sum K''K'' = 2 e^-2y (I0 - I1)(2y)', sKpKp, &
                        2._xp*EXP(-2._xp*y)*(bessel_i(0,2._xp*y)-bessel_i(1,2._xp*y)), tol)
      ENDDO
    ENDDO
  ELSE
    CALL section('the coupled operator is never built when BPAR is off')
    CALL check('bp_phi_rho stays zero', ALL(bp_phi_rho .EQ. 0._xp))
    CALL check('bp_xcoupl stays zero',  ALL(bp_xcoupl  .EQ. 0._xp))
    CALL check('bp_dB_pp stays zero',   ALL(bp_dB_pp   .EQ. 0._xp))
  ENDIF

  !------------------------------------------------------ field solve
  CALL section('the fields solved from random moments satisfy the field equations')
  CALL random_p0_moments
  CALL solve_EM_fields
  CALL compute_nadiab_moments
  ! Residuals of quasineutrality and of perpendicular pressure balance, from
  ! the non adiabatic moments h the code built. With BPAR off the second one is
  ! simply the statement that dBpar vanishes.
  worst = 0._xp; scale = 0._xp
  DO ikx = 1,local_nkx
    DO iky = 1,local_nky
      IF ((ikx .EQ. ikx0) .AND. (iky .EQ. iky0)) CYCLE
      R_qn = 0._xp; R_B = hatB(izi,ieven)**2*dBpar(iky,ikx,izi)
      DO ia = 1,local_na
        rho = 0._xp; pperp = 0._xp
        DO ij = 1,local_nj
          iji = ij + ngj_o2
          hj = nadiab_p0(ia,iji,iky,ikx)
          rho   = rho   + kernel   (ia,iji,iky,ikx,izi,ieven)*hj
          pperp = pperp + kernel_dB(ia,iji,iky,ikx,izi,ieven)*hj
        ENDDO
        R_qn = R_qn + q(ia)*(rho - q(ia)/tau(ia)*phi(iky,ikx,izi))
        IF (BPAR) R_B = R_B + 0.5_xp*beta*tau(ia)*pperp
      ENDDO
      worst = MAX(worst, ABS(R_qn), ABS(R_B))
      scale = MAX(scale, ABS(phi(iky,ikx,izi)), ABS(dBpar(iky,ikx,izi)))
    ENDDO
  ENDDO
  CALL check('the fields are not identically zero', scale .GT. 0._xp)
  CALL check('max residual of both field equations is round off', worst .LT. tol*scale)
  CALL check('phi is zero at the origin',   phi(iky0,ikx0,izi)   .EQ. (0._xp,0._xp))
  CALL check('dBpar is zero at the origin', dBpar(iky0,ikx0,izi) .EQ. (0._xp,0._xp))
  IF (.NOT. BPAR) THEN
    CALL section('with BPAR off the A_par-only path is untouched')
    CALL check('dBpar is identically zero', ALL(dBpar .EQ. (0._xp,0._xp)))
    worst = 0._xp
    DO ikx = 1,local_nkx
      DO iky = 1,local_nky
        rho = 0._xp
        DO ij = 1,local_nj
          iji = ij + ngj_o2
          DO ia = 1,local_na
            rho = rho + q(ia)*kernel(ia,iji,iky,ikx,izi,ieven)*moments(ia,ip0,iji,iky,ikx,izi,updatetlevel)
          ENDDO
        ENDDO
        worst = MAX(worst, ABS(phi(iky,ikx,izi) - inv_poisson_op(iky,ikx,1)*rho))
      ENDDO
    ENDDO
    CALL check('phi is the electrostatic Poisson solution', worst .LT. tol*scale)
  ELSE
    CALL section('dBpar is a correction of order beta')
    ! Pressure balance bounds |dBpar| by (beta/2)(|pperp| + |C phi|)/M with
    ! M >= 1; a looser but convention free statement is enough here.
    worst = 0._xp
    DO ikx = 1,local_nkx
      DO iky = 1,local_nky
        worst = MAX(worst, ABS(dBpar(iky,ikx,izi)))
      ENDDO
    ENDDO
    CALL check('max |dBpar| < beta * max(|phi|,|dBpar|)', worst .LT. beta*scale)
  ENDIF

  !------------------------------------------------------ gradient drive
  CALL section('the hierarchy drive of dBpar alone')
  ! With k_cB = k_gB = 0 (set in the deck), a single z plane, no collisions and
  ! no hyperdiffusion, the RHS of a state with only dBpar nonzero is exactly
  ! -DdB = -i ky (tau/q)(xphij K'_j + xphijp1 K'_{j+1} + xphijm1 K'_{j-1}) dBpar
  ! for p = 0 and 2, and zero for every other p. With BPAR off the field has to
  ! be ignored altogether.
  moments = 0._xp; phi = 0._xp; psi = 0._xp
  dBpar = B_AMP
  CALL compute_nadiab_moments
  CALL compute_gradients_z
  CALL compute_interp_z
  CALL compute_Sapj
  CALL compute_moments_eq_rhs
  worst = 0._xp
  DO ikx = 1,local_nkx
    DO iky = 1,local_nky
      DO ia = 1,local_na
        y = sigma2_tau_o2(ia)*kp2array(iky,ikx,izi,ieven)
        CALL closed_form_kernels(y, Kc)
        DO ij = 1,local_nj
          j_int = jarray(ij+ngj_o2)
          DO ip = 1,local_np
            p_int = parray(ip+ngp_o2)
            got  = moments_rhs(ia,ip,ij,iky,ikx,1,updatetlevel)
            want = 0._xp
            IF (BPAR) THEN
              IF (p_int .EQ. 0) THEN
                want = -imagu*kyarray(iky)*tau(ia)/q(ia)*B_AMP*( &
                       (k_N(ia) + 2._xp*j_int*k_T(ia))*kprime(Kc, j_int) &
                      - k_T(ia)*(j_int+1._xp)         *kprime(Kc, j_int+1) &
                      - k_T(ia)* j_int                *kprime(Kc, j_int-1))
              ELSEIF (p_int .EQ. 2) THEN
                want = -imagu*kyarray(iky)*tau(ia)/q(ia)*B_AMP*k_T(ia)/SQRT2*kprime(Kc, j_int)
              ENDIF
            ENDIF
            worst = MAX(worst, ABS(got - want))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  CALL check('max deviation of the RHS from the closed form drive is round off', worst .LT. tol)
  IF (BPAR) THEN
    ! and the drift kinetic limit is the pressure gradient drive (k_N + k_T)
    iky = 2; ikx = ikx0; ia = 1
    y = sigma2_tau_o2(ia)*kp2array(iky,ikx,izi,ieven)
    got  = moments_rhs(ia,1+(ip0-1-ngp_o2),ij0-ngj_o2,iky,ikx,1,updatetlevel)
    want = -imagu*kyarray(iky)*tau_q(ia)*B_AMP*EXP(-y)*(k_N(ia) + k_T(ia)*(1._xp - y))
    CALL check_cmpx('p = 0, j = 0 drive is -i ky (tau/q) e^-y (k_N + k_T (1-y)) dBpar', got, want, tol)
  ENDIF

  CALL finish_suite

CONTAINS

  !> K_j(y) = y^j e^-y / j! for j = 0..jmax+1, with K_{-1} = 0.
  SUBROUTINE closed_form_kernels(y_, K_)
    REAL(xp), INTENT(IN)  :: y_
    REAL(xp), INTENT(OUT) :: K_(-1:T_JMAX+1)
    INTEGER :: j_
    K_(-1) = 0._xp
    K_(0)  = EXP(-y_)
    DO j_ = 1,T_JMAX+1
      K_(j_) = K_(j_-1)*y_/REAL(j_,xp)
    ENDDO
  END SUBROUTINE closed_form_kernels

  !> K'_j = K_j - K_{j-1} of the closed form kernels, zero below j = 0 like the
  !! ghost of kernel_dB.
  REAL(xp) FUNCTION kprime(K_, j_)
    REAL(xp), INTENT(IN) :: K_(-1:T_JMAX+1)
    INTEGER,  INTENT(IN) :: j_
    IF (j_ .LT. 0) THEN
      kprime = 0._xp
    ELSE
      kprime = K_(j_) - K_(j_-1)
    ENDIF
  END FUNCTION kprime

  !> Laguerre polynomials L_0..L_{jmax+1} at x by the three term recurrence.
  SUBROUTINE laguerre(x_, L_)
    REAL(xp), INTENT(IN)  :: x_
    REAL(xp), INTENT(OUT) :: L_(0:T_JMAX+1)
    INTEGER :: n_
    L_(0) = 1._xp
    L_(1) = 1._xp - x_
    DO n_ = 1,T_JMAX
      L_(n_+1) = ((2._xp*n_ + 1._xp - x_)*L_(n_) - n_*L_(n_-1))/REAL(n_+1,xp)
    ENDDO
  END SUBROUTINE laguerre

  !> Modified Bessel function I_n(z) from its power series (z <= 4 here).
  REAL(xp) FUNCTION bessel_i(n_, z_)
    INTEGER,  INTENT(IN) :: n_
    REAL(xp), INTENT(IN) :: z_
    REAL(xp) :: term
    INTEGER  :: k_
    term = (0.5_xp*z_)**n_/GAMMA(REAL(n_+1,xp))
    bessel_i = term
    DO k_ = 1,60
      term = term*(0.5_xp*z_)**2/(REAL(k_,xp)*REAL(k_+n_,xp))
      bessel_i = bessel_i + term
      IF (term .LT. EPSILON(1._xp)*bessel_i) EXIT
    ENDDO
  END FUNCTION bessel_i

  !> The p = 0 non adiabatic moment the code built, read back from its array.
  COMPLEX(xp) FUNCTION nadiab_p0(ia_, iji_, iky_, ikx_)
    USE array, ONLY: nadiab_moments
    INTEGER, INTENT(IN) :: ia_, iji_, iky_, ikx_
    nadiab_p0 = nadiab_moments(ia_,ip0,iji_,iky_,ikx_,izi)
  END FUNCTION nadiab_p0

  !> Random p = 0 moments on the interior j, zero elsewhere (ghosts, p > 0).
  SUBROUTINE random_p0_moments
    REAL(xp) :: re, im
    INTEGER  :: ia_, ij_, iky_, ikx_
    moments = 0._xp
    DO ikx_ = 1,local_nkx
      DO iky_ = 1,local_nky
        DO ij_ = 1,local_nj
          DO ia_ = 1,local_na
            CALL RANDOM_NUMBER(re); CALL RANDOM_NUMBER(im)
            moments(ia_,ip0,ij_+ngj_o2,iky_,ikx_,izi,updatetlevel) = CMPLX(re-0.5_xp, im-0.5_xp, xp)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE random_p0_moments

  !> Feed the namelists through a scratch file so nothing is left on disk.
  SUBROUTINE write_input_deck
    OPEN(UNIT=lu_in, STATUS='SCRATCH', ACTION='READWRITE', FORM='FORMATTED')
    WRITE(lu_in,'(A)')       '&GRID'
    WRITE(lu_in,'(A,I0)')    '  pmax = ', T_PMAX
    WRITE(lu_in,'(A,I0)')    '  jmax = ', T_JMAX
    WRITE(lu_in,'(A,I0)')    '  Nx   = ', T_NX
    WRITE(lu_in,'(A,F8.2)')  '  Lx   = ', T_LX
    WRITE(lu_in,'(A,I0)')    '  Ny   = ', T_NY
    WRITE(lu_in,'(A,F8.2)')  '  Ly   = ', T_LY
    WRITE(lu_in,'(A,I0)')    '  Nz   = ', T_NZ
    WRITE(lu_in,'(A)')       '  SG   = .false.'
    WRITE(lu_in,'(A)')       '  Nexc = 1'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&GEOMETRY'
    WRITE(lu_in,'(A)')       "  geom = 'z-pinch'"
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&DIAGNOSTICS'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&MODEL'
    WRITE(lu_in,'(A)')       "  LINEARITY = 'linear'"
    WRITE(lu_in,'(A)')       '  Na      = 2'
    WRITE(lu_in,'(A)')       '  ADIAB_E = .false.'
    WRITE(lu_in,'(A,ES23.15)')'  beta    = ', beta_in
    WRITE(lu_in,'(A,L1)')    '  BPAR    = ', bpar_wanted
    WRITE(lu_in,'(A)')       '  nu      = 0.0'
    ! no curvature and no grad B: the perpendicular magnetic term of the
    ! hierarchy vanishes and the RHS of the drive check is the drive alone
    WRITE(lu_in,'(A)')       '  k_cB    = 0.0'
    WRITE(lu_in,'(A)')       '  k_gB    = 0.0'
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
    WRITE(lu_in,'(A)')       '  k_N_   = 2.0'
    WRITE(lu_in,'(A)')       '  k_T_   = 0.4'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&SPECIES'
    WRITE(lu_in,'(A)')       "  name_  = 'electrons'"
    WRITE(lu_in,'(A)')       '  tau_   = 0.8'
    WRITE(lu_in,'(A)')       '  sigma_ = 0.2'
    WRITE(lu_in,'(A)')       '  q_     = -1.0'
    WRITE(lu_in,'(A)')       '  k_N_   = 2.0'
    WRITE(lu_in,'(A)')       '  k_T_   = 1.5'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&COLLISION'
    WRITE(lu_in,'(A)')       "  collision_model = 'DG'"
    WRITE(lu_in,'(A)')       '  GK_CO           = .true.'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&INITIAL'
    WRITE(lu_in,'(A)')       "  INIT_OPT = 'phi'"
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&TIME_INTEGRATION'
    WRITE(lu_in,'(A)')       "  numerical_scheme = 'RK4'"
    WRITE(lu_in,'(A)')       '/'
    REWIND(lu_in)
  END SUBROUTINE write_input_deck

END PROGRAM test_bpar
