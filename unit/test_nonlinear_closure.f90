PROGRAM test_nonlinear_closure
  ! Exercises the Laguerre truncation of the nonlinear term (closure_mod,
  ! nonlinear_mod): which couplings (n,j,s) of
  !   S_apj = sum_n { K_n phi , sum_s d_njs N_aps }
  ! each nonlinear_closure retains, and whether the resulting ExB and flutter
  ! brackets conserve the free energy sum_pj |N_apj|^2 of every species.
  !
  ! The conservation check is the one that matters in practice: the dealiased
  ! bracket is antisymmetric, so a coupling (n,j,s) and its mirror (n,s,j)
  ! exchange energy without creating any. A closure that keeps a coupling
  ! without its mirror pumps energy into N_apj. On random fields this shows up
  ! as an O(1) fraction of sum |N||S|, while a symmetric closure gives round off.
  !
  ! The closure is taken from the command line because readinputs can only
  ! run once per process (allocate_array does not deallocate), so "make check"
  ! runs this binary once per closure:
  !
  !   full_sum                    conserving
  !   anti_laguerre_aliasing      not conserving
  !   anti_laguerre_aliasing_sym  conserving
  !   truncation                  nmax = -1, i.e. anti_laguerre_aliasing, not conserving
  !   truncation_nmax1            nmax =  1, uniform bound, conserving
  !
  ! Single rank, z-pinch with a single z plane and two kinetic species, beta > 0
  ! so that the flutter bracket is exercised too.
  USE prec_const, ONLY: xp
  USE basic,      ONLY: lu_in
  USE parallel,   ONLY: ppinit
  USE grid,       ONLY: total_nkx, total_nky, local_nkx, local_nky, local_na, local_np, local_nj, &
                        ngj_o2, ngp_o2, ngz_o2, jarray, kxarray, kx_max, iky0, jmax, pmax
  USE fields,     ONLY: moments, phi, psi
  USE array,      ONLY: Sapj, dnjs
  USE model,      ONLY: EM
  USE closure,    ONLY: nonlinear_closure, nmax, nmaxarray, smaxarray
  USE time_integration, ONLY: set_updatetlevel
  USE nonlinear,  ONLY: compute_Sapj, nonlinear_init
  USE test_harness
  IMPLICIT NONE

  INTEGER,  PARAMETER :: T_PMAX = 3, T_JMAX = 3
  INTEGER,  PARAMETER :: T_NX = 16, T_NY = 16, T_NZ = 1
  REAL(xp), PARAMETER :: T_LX = 40._xp, T_LY = 40._xp
  ! Conserving closures land at ~1e-17 of sum |N||S|; the unpaired couplings of the
  ! others leak 1e-3 to 5e-2 of it on the fields below (the flutter terms, which
  ! only exchange energy between Hermite degrees, inflate the scale of the total).
  REAL(xp), PARAMETER :: NONCONS_FLOOR  = 1.0e-4_xp
  REAL(xp), PARAMETER :: EXCHANGE_FLOOR = 1.0e-3_xp ! per-p exchange by the flutter bracket

  CHARACTER(len=32) :: mode, closure_in
  INTEGER  :: nmax_in, nargs, ia, ip, ij, in, is, n, j, s, iz
  LOGICAL  :: conserving, sym_ok, rule_ok, in_range, coupled
  REAL(xp) :: tol, inj, scale, worst_rel, Pinj(0:T_PMAX), Pscale(0:T_PMAX)
  LOGICAL, ALLOCATABLE :: kept(:,:,:)

  CALL ppinit
  nargs = COMMAND_ARGUMENT_COUNT()
  IF (nargs .GE. 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, mode)
  ELSE
    mode = 'full_sum'
  ENDIF
  SELECT CASE (TRIM(mode))
  CASE ('full_sum');                   closure_in = 'full_sum';                   nmax_in = -1; conserving = .TRUE.
  CASE ('anti_laguerre_aliasing');     closure_in = 'anti_laguerre_aliasing';     nmax_in = -1; conserving = .FALSE.
  CASE ('anti_laguerre_aliasing_sym'); closure_in = 'anti_laguerre_aliasing_sym'; nmax_in = -1; conserving = .TRUE.
  CASE ('truncation');                 closure_in = 'truncation';                 nmax_in = -1; conserving = .FALSE.
  CASE ('truncation_nmax1');           closure_in = 'truncation';                 nmax_in =  1; conserving = .TRUE.
  CASE DEFAULT
    WRITE(*,'(A)') 'unknown mode "'//TRIM(mode)//'"'
    CALL start_suite('nonlinear closure ['//TRIM(mode)//']')
    CALL check('mode is known to the test', .FALSE.)
    CALL finish_suite
    STOP
  END SELECT

  CALL start_suite('nonlinear closure ['//TRIM(mode)//']')
  tol = tol_of(1.0e5_xp)
  iz  = 1

  CALL write_input_deck
  CALL readinputs
  CALL auxval
  CALL nonlinear_init ! work arrays of the bracket, allocated by control in a run
  CALL set_updatetlevel(1)

  !--------------------------------------------------------------- Setup
  CALL section('the run was configured the way the test assumes')
  CALL check_int('single rank: local_nky is the total', local_nky, total_nky)
  CALL check_int('single rank: local_nkx is the total', local_nkx, total_nkx)
  CALL check_int('two kinetic species', local_na, 2)
  CALL check_int('pmax as requested', pmax, T_PMAX)
  CALL check_int('jmax as requested', jmax, T_JMAX)
  CALL check_int('all Laguerre degrees on this rank', local_nj, T_JMAX+1)
  CALL check('closure name read back', TRIM(nonlinear_closure) .EQ. TRIM(closure_in))
  CALL check_int('nmax read back', nmax, nmax_in)
  CALL check('electromagnetic (flutter bracket active)', EM)

  !----------------------------------------------------- retained couplings
  ! kept(n,j,s): the coupling N_pj <- {K_n phi, N_ps} is in the sum, for the
  ! triplets the Laguerre triple product does not kill anyway.
  CALL section('retained couplings follow the rule of the closure')
  ALLOCATE(kept(0:T_JMAX,0:T_JMAX,0:T_JMAX))
  kept = .FALSE.
  DO ij = 1,local_nj
    j = jarray(ij+ngj_o2)
    DO in = 1,T_JMAX+1
      n = in-1
      IF (n .GT. nmaxarray(ij)) CYCLE
      DO is = 1,T_JMAX+1
        s = is-1
        IF (s .GT. smaxarray(in,ij)) CYCLE
        kept(n,j,s) = ABS(dnjs(in,ij,is)) .GT. 0._xp
      ENDDO
    ENDDO
  ENDDO
  rule_ok = .TRUE.
  DO n = 0,T_JMAX
    DO j = 0,T_JMAX
      DO s = 0,T_JMAX
        coupled = ABS(dnjs(n+1,j+1,s+1)) .GT. 0._xp
        SELECT CASE (TRIM(mode))
        CASE ('full_sum')
          rule_ok = rule_ok .AND. (kept(n,j,s) .EQV. coupled)
        CASE ('anti_laguerre_aliasing', 'truncation')
          rule_ok = rule_ok .AND. (kept(n,j,s) .EQV. (coupled .AND. (n+j .LE. T_JMAX)))
        CASE ('anti_laguerre_aliasing_sym')
          rule_ok = rule_ok .AND. (kept(n,j,s) .EQV. (coupled .AND. (n+MAX(j,s) .LE. T_JMAX)))
        CASE ('truncation_nmax1')
          rule_ok = rule_ok .AND. (kept(n,j,s) .EQV. (coupled .AND. (n .LE. nmax_in)))
        END SELECT
      ENDDO
    ENDDO
  ENDDO
  CALL check('kept(n,j,s) matches the documented rule', rule_ok)
  ! the triple product itself: the Laguerre triangle rule, which the symmetry
  ! argument of the conservation check relies on
  rule_ok = .TRUE.; sym_ok = .TRUE.; in_range = .TRUE.
  DO n = 0,T_JMAX
    DO j = 0,T_JMAX
      DO s = 0,T_JMAX
        rule_ok  = rule_ok  .AND. (ABS(dnjs(n+1,j+1,s+1) - dnjs(n+1,s+1,j+1)) .EQ. 0._xp)
        sym_ok   = sym_ok   .AND. (kept(n,j,s) .EQV. kept(n,s,j))
        in_range = in_range .AND. ((.NOT. kept(n,j,s)) .OR. (n + MAX(j,s) .LE. T_JMAX))
      ENDDO
    ENDDO
  ENDDO
  CALL check('d_njs is symmetric in (j,s)', rule_ok)
  IF (conserving) THEN
    CALL check('the set of retained couplings is symmetric in (j,s)', sym_ok)
  ELSE
    CALL check('the set of retained couplings is NOT symmetric in (j,s)', .NOT. sym_ok)
    ! and name one unpaired coupling so a reader sees what is at stake
    CALL check('N_p0 <- {K_J phi, N_pJ} is kept',      kept(T_JMAX,0,T_JMAX))
    CALL check('its mirror N_pJ <- {K_J phi, N_p0} is dropped', .NOT. kept(T_JMAX,T_JMAX,0))
  ENDIF
  IF (TRIM(mode) .EQ. 'anti_laguerre_aliasing_sym') THEN
    CALL check('no retained product L_n L_j or L_n L_s leaves the Laguerre range', in_range)
  ENDIF

  !---------------------------------------------- free energy, ES bracket only
  ! With psi = 0 the bracket does not mix Hermite degrees, so the invariant is
  ! sum_j |N_pj|^2 for every (species, p) separately.
  CALL section('ExB bracket: sum over j and k of Re(N* S) per species and p')
  CALL fill_random_fields(with_psi=.FALSE.)
  CALL compute_Sapj
  worst_rel = 0._xp
  DO ia = 1,local_na
    CALL injection_per_p(ia, Pinj, Pscale)
    DO ip = 0,T_PMAX
      worst_rel = MAX(worst_rel, ABS(Pinj(ip))/Pscale(ip))
    ENDDO
    WRITE(*,'(A,I0,A,*(ES10.2))') '  species ', ia, ' |sum Re(N*S)|/sum|N||S| per p: ', ABS(Pinj)/Pscale
  ENDDO
  IF (conserving) THEN
    CALL check('every (a,p) sum vanishes to round off', worst_rel .LT. tol)
  ELSE
    CALL check('some (a,p) sum is a finite fraction of sum |N||S|', worst_rel .GT. NONCONS_FLOOR)
  ENDIF
  CALL check('the bracket is not trivially zero', MAXVAL(ABS(Sapj)) .GT. 0._xp)

  !------------------------------------- free energy, ES + flutter brackets
  ! The flutter bracket couples p to p+-1 with the same Laguerre sums; the
  ! invariant is then the full sum_pj |N_pj|^2 of each species.
  CALL section('ExB + flutter brackets: sum over p, j and k of Re(N* S) per species')
  CALL fill_random_fields(with_psi=.TRUE.)
  CALL compute_Sapj
  worst_rel = 0._xp
  DO ia = 1,local_na
    CALL injection_per_p(ia, Pinj, Pscale)
    inj = SUM(Pinj); scale = SUM(Pscale)
    worst_rel = MAX(worst_rel, ABS(inj)/scale)
    WRITE(*,'(A,I0,A,ES10.2,A,*(ES10.2))') '  species ', ia, ' |sum Re(N*S)|/sum|N||S| total: ', ABS(inj)/scale, '  per p: ', ABS(Pinj)/Pscale
    ! the flutter part alone redistributes energy between p's: with psi on the
    ! per-p sums must not vanish even for a conserving closure
    IF (ia .EQ. 1 .AND. conserving) &
      CALL check('with psi the per-p sums are not individually zero', &
                 MAXVAL(ABS(Pinj)/Pscale) .GT. EXCHANGE_FLOOR)
  ENDDO
  IF (conserving) THEN
    CALL check('every species sum vanishes to round off', worst_rel .LT. tol)
  ELSE
    CALL check('some species sum is a finite fraction of sum |N||S|', worst_rel .GT. NONCONS_FLOOR)
  ENDIF

  CALL finish_suite

CONTAINS

  !> Random moments, phi and (optionally) psi on the evolved (p,j) range, with
  !! the ky = 0 row made Hermitian in kx so that the half spectrum describes a
  !! real field (the c2r transforms of the bracket assume it).
  SUBROUTINE fill_random_fields(with_psi)
    LOGICAL, INTENT(IN) :: with_psi
    INTEGER :: ia, ip, ij, iky, ikx, izi, ipi, iji, nseed
    INTEGER, ALLOCATABLE :: seed(:)
    REAL(xp) :: r(2)
    CALL RANDOM_SEED(SIZE=nseed)
    ALLOCATE(seed(nseed)); seed = 20261027
    CALL RANDOM_SEED(PUT=seed)
    izi = iz + ngz_o2
    moments = 0._xp; phi = 0._xp; psi = 0._xp
    DO ia = 1,local_na
      DO ip = 1,local_np
        ipi = ip + ngp_o2
        DO ij = 1,local_nj
          iji = ij + ngj_o2
          DO ikx = 1,local_nkx
            DO iky = 1,local_nky
              CALL RANDOM_NUMBER(r)
              moments(ia,ipi,iji,iky,ikx,izi,1) = CMPLX(r(1)-0.5_xp, r(2)-0.5_xp, xp)
            ENDDO
          ENDDO
          CALL hermitian_ky0_row(moments(ia,ipi,iji,:,:,izi,1))
        ENDDO
      ENDDO
    ENDDO
    DO ikx = 1,local_nkx
      DO iky = 1,local_nky
        CALL RANDOM_NUMBER(r)
        phi(iky,ikx,izi) = CMPLX(r(1)-0.5_xp, r(2)-0.5_xp, xp)
        IF (with_psi) THEN
          CALL RANDOM_NUMBER(r)
          psi(iky,ikx,izi) = 0.3_xp*CMPLX(r(1)-0.5_xp, r(2)-0.5_xp, xp)
        ENDIF
      ENDDO
    ENDDO
    CALL hermitian_ky0_row(phi(:,:,izi))
    IF (with_psi) CALL hermitian_ky0_row(psi(:,:,izi))
  END SUBROUTINE fill_random_fields

  !> F(kx, ky=0) = conjg(F(-kx, ky=0)), F(0,0) real, modes without a partner zeroed.
  SUBROUTINE hermitian_ky0_row(F)
    COMPLEX(xp), INTENT(INOUT) :: F(:,:) ! (ky, kx)
    INTEGER  :: ikx, jkx
    REAL(xp) :: kx
    DO ikx = 1,local_nkx
      kx = kxarray(iky0,ikx)
      IF (kx .EQ. 0._xp) THEN
        F(iky0,ikx) = CMPLX(REAL(F(iky0,ikx),xp), 0._xp, xp)
      ELSE IF (kx .GT. 0._xp) THEN
        jkx = 0
        DO is = 1,local_nkx
          IF (ABS(kxarray(iky0,is) + kx) .LT. 1.0e-10_xp*kx_max) jkx = is
        ENDDO
        IF (jkx .EQ. 0) THEN
          F(iky0,ikx) = 0._xp
        ELSE
          F(iky0,jkx) = CONJG(F(iky0,ikx))
        ENDIF
      ENDIF
    ENDDO
    ! the negative Nyquist mode has no positive partner
    DO ikx = 1,local_nkx
      kx = kxarray(iky0,ikx)
      IF (kx .LT. 0._xp) THEN
        jkx = 0
        DO is = 1,local_nkx
          IF (ABS(kxarray(iky0,is) + kx) .LT. 1.0e-10_xp*kx_max) jkx = is
        ENDDO
        IF (jkx .EQ. 0) F(iky0,ikx) = 0._xp
      ENDIF
    ENDDO
  END SUBROUTINE hermitian_ky0_row

  !> For species ia: per Hermite degree, the half-spectrum sums
  !!   inj(p)   = sum_{j,k} w_k Re( conjg(N_pj) S_pj ),   scale(p) = sum_{j,k} w_k |N_pj| |S_pj|,
  !! with w = 1 on the ky = 0 row and 2 elsewhere (the conjugate half).
  SUBROUTINE injection_per_p(ia, inj, scale)
    INTEGER,  INTENT(IN)  :: ia
    REAL(xp), INTENT(OUT) :: inj(0:T_PMAX), scale(0:T_PMAX)
    INTEGER  :: ip, ij, iky, ikx, izi, ipi, iji
    REAL(xp) :: w
    COMPLEX(xp) :: Nk, Sk
    izi = iz + ngz_o2
    inj = 0._xp; scale = 0._xp
    DO ip = 1,local_np
      ipi = ip + ngp_o2
      DO ij = 1,local_nj
        iji = ij + ngj_o2
        DO ikx = 1,local_nkx
          DO iky = 1,local_nky
            w = 2._xp
            IF (iky .EQ. iky0) w = 1._xp
            Nk = moments(ia,ipi,iji,iky,ikx,izi,1)
            Sk = Sapj(ia,ip,ij,iky,ikx,iz)
            inj(ip-1)   = inj(ip-1)   + w*REAL(CONJG(Nk)*Sk, xp)
            scale(ip-1) = scale(ip-1) + w*ABS(Nk)*ABS(Sk)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE injection_per_p

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
    WRITE(lu_in,'(A)')       "  LINEARITY = 'nonlinear'"
    WRITE(lu_in,'(A)')       '  Na      = 2'
    WRITE(lu_in,'(A)')       '  ADIAB_E = .false.'
    WRITE(lu_in,'(A)')       '  beta    = 0.05'
    WRITE(lu_in,'(A)')       '  nu      = 0.0'
    WRITE(lu_in,'(A)')       '  k_cB    = 0.0'
    WRITE(lu_in,'(A)')       '  k_gB    = 0.0'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&CLOSURE'
    WRITE(lu_in,'(A)')       "  hierarchy_closure = 'truncation'"
    WRITE(lu_in,'(A)')       "  nonlinear_closure = '"//TRIM(closure_in)//"'"
    WRITE(lu_in,'(A,I0)')    '  nmax = ', nmax_in
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

END PROGRAM test_nonlinear_closure
