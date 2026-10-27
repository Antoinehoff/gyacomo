PROGRAM test_adiab_e
  ! Exercises the adiabatic electron response of the quasineutrality equation,
  !   rho = P phi + (phi - <phi>),
  ! with P(z) the ion polarisation and <.> the flux surface average, which the
  ! ky = 0 modes close with the exact form
  !   <phi> = <phi0>/(1 - <1/(P+1)>),   phi0 = rho/(P+1),
  ! and the ky /= 0 modes with <phi> = 0 (solve_EM_fields::poisson,
  ! numerics_mod::evaluate_poisson_op).
  !
  ! The polarisation and the densities are rebuilt here from the code's kernels
  ! and moments, so the test pins the solve to what the hierarchy is consistent
  ! with, whatever the kernel model; the average is rebuilt from the Simpson
  ! weights and the Jacobian rather than taken from the code's helpers.
  !
  ! Single rank, s-alpha with a large inverse aspect ratio so that hatB, hence
  ! P(z), varies strongly along the field line: this is where the former
  ! <phi> = <rho/P> shortcut (exact only for a z-independent P) goes wrong.
  USE prec_const, ONLY: xp
  USE basic,      ONLY: lu_in
  USE parallel,   ONLY: ppinit
  USE grid,       ONLY: total_nkx, total_nky, local_nkx, local_nky, local_na, local_nz, local_nj, &
                        ngj_o2, ngz_o2, ikx0, iky0, ip0, kyarray, kxarray, zweights_SR, ieven
  USE fields,     ONLY: moments, phi
  USE array,      ONLY: kernel, inv_poisson_op, inv_poisson_op_fsa
  USE model,      ONLY: ADIAB_E
  USE species,    ONLY: q, q2_tau
  USE geometry,   ONLY: hatB, Jacobian
  USE time_integration, ONLY: updatetlevel, set_updatetlevel
  USE test_harness
  IMPLICIT NONE

  INTEGER,  PARAMETER :: T_PMAX = 2, T_JMAX = 4
  INTEGER,  PARAMETER :: T_NX = 8, T_NY = 4, T_NZ = 16
  REAL(xp), PARAMETER :: T_LX = 40._xp, T_LY = 40._xp
  REAL(xp), PARAMETER :: T_EPS = 0.3_xp

  REAL(xp) :: tol, worst, worst_e, worst_c, worst_k, scale, Bmin, Bmax, sw, pz
  REAL(xp), DIMENSION(T_NZ) :: P, w
  COMPLEX(xp), DIMENSION(T_NZ) :: rho, r
  COMPLEX(xp) :: c, fsa_phi, fsa_e, fsa_num, fsa_old, got
  INTEGER  :: ikx, iky, iz, izi, ij, iji, ia
  LOGICAL  :: used_fsa

  CALL ppinit
  CALL start_suite('adiabatic electron response')
  tol = tol_of(1.0e4_xp)

  CALL write_input_deck
  CALL readinputs
  CALL auxval
  CALL set_updatetlevel(1) ! the time level the solver reads, set by initial_mod in a run

  !--------------------------------------------------------------- Setup
  CALL section('the run was configured the way the test assumes')
  CALL check_int('single rank: local_nky is the total', local_nky, total_nky)
  CALL check_int('single rank: local_nkx is the total', local_nkx, total_nkx)
  CALL check_int('single rank: local_nz is Nz', local_nz, T_NZ)
  CALL check_int('one kinetic species', local_na, 1)
  CALL check('adiabatic electrons are on', ADIAB_E)
  CALL check('the origin is on this rank', (ikx0 .GE. 1) .AND. (iky0 .GE. 1))
  Bmin = HUGE(1._xp); Bmax = 0._xp; sw = 0._xp
  DO iz = 1,local_nz
    izi = iz + ngz_o2
    Bmin = MIN(Bmin, hatB(izi,ieven)); Bmax = MAX(Bmax, hatB(izi,ieven))
    w(iz) = zweights_SR(iz)*Jacobian(izi,ieven)
    sw = sw + w(iz)
  ENDDO
  w = w/sw ! normalised flux surface average weights
  CALL check('hatB varies by more than 50 % along the field line', Bmax/Bmin .GT. 1.5_xp)
  CALL check('the average weights are positive', ALL(w .GT. 0._xp))

  !------------------------------------------------- precomputed <1/(P+1)>
  CALL section('evaluate_poisson_op stores <1/(P+1)> of the ky = 0 row')
  worst = 0._xp
  iky = iky0
  DO ikx = 1,local_nkx
    IF (ikx .EQ. ikx0) CYCLE
    CALL polarisation(iky, ikx, P)
    got = 0._xp
    DO iz = 1,local_nz
      got = got + w(iz)/(P(iz) + 1._xp)
    ENDDO
    worst = MAX(worst, ABS(REAL(got,xp) - inv_poisson_op_fsa(ikx)))
    ! and the pointwise operator is 1/(P+1) for kx /= 0
    DO iz = 1,local_nz
      worst = MAX(worst, ABS(inv_poisson_op(iky,ikx,iz)*(P(iz) + 1._xp) - 1._xp))
    ENDDO
  ENDDO
  CALL check('max deviation from the rebuilt <1/(P+1)> and 1/(P+1) is round off', worst .LT. tol)
  CALL check('it is zero at the origin', inv_poisson_op_fsa(ikx0) .EQ. 0._xp)
  CALL check('it lies strictly inside (0, 1) away from the origin', &
             ALL(inv_poisson_op_fsa(1:ikx0-1) .GT. 0._xp .AND. inv_poisson_op_fsa(1:ikx0-1) .LT. 1._xp) .AND. &
             ALL(inv_poisson_op_fsa(ikx0+1:local_nkx) .GT. 0._xp .AND. inv_poisson_op_fsa(ikx0+1:local_nkx) .LT. 1._xp))

  !------------------------------------------------------ field solve, ky = 0
  CALL section('the ky = 0 solve from random moments satisfies rho = P phi + phi - <phi>')
  CALL random_p0_moments
  CALL solve_EM_fields
  worst = 0._xp; worst_e = 0._xp; worst_c = 0._xp; worst_k = 0._xp; scale = 0._xp
  used_fsa = .FALSE.
  iky = iky0
  DO ikx = 1,local_nkx
    IF (ikx .EQ. ikx0) CYCLE
    CALL polarisation(iky, ikx, P)
    CALL charge_density(iky, ikx, rho)
    ! r(z) = rho - (P+1) phi has to be the z-independent constant -<phi>
    fsa_phi = 0._xp; fsa_e = 0._xp; fsa_num = 0._xp; fsa_old = 0._xp; pz = 0._xp
    DO iz = 1,local_nz
      izi = iz + ngz_o2
      r(iz)   = rho(iz) - (P(iz) + 1._xp)*phi(iky,ikx,izi)
      fsa_phi = fsa_phi + w(iz)*phi(iky,ikx,izi)                 ! <phi>
      fsa_e   = fsa_e   + w(iz)*(rho(iz) - P(iz)*phi(iky,ikx,izi)) ! <n_e> = <phi - <phi>>
      fsa_num = fsa_num + w(iz)*rho(iz)/(P(iz) + 1._xp)          ! <rho/(P+1)>
      pz      = pz      + w(iz)*P(iz)/(P(iz) + 1._xp)            ! <P/(P+1)>
      fsa_old = fsa_old + w(iz)*rho(iz)/P(iz)                    ! the former <rho/P>
      scale   = MAX(scale, ABS(phi(iky,ikx,izi)))
    ENDDO
    c = -SUM(r)/REAL(local_nz,xp)
    DO iz = 1,local_nz
      worst = MAX(worst, ABS(r(iz) + c))
    ENDDO
    worst_e = MAX(worst_e, ABS(fsa_e))
    worst_c = MAX(worst_c, ABS(c - fsa_phi), ABS(c - fsa_num/pz))
    worst_k = MAX(worst_k, ABS(fsa_old - c))
    used_fsa = used_fsa .OR. (ABS(c) .GT. 0._xp)
  ENDDO
  CALL check('the fields are not identically zero', scale .GT. 0._xp)
  CALL check('rho - (P+1) phi is constant along the field line (round off)', worst .LT. tol*scale)
  CALL check('that constant is -<phi>, the flux surface average with the Jacobian', worst_c .LT. tol*scale)
  CALL check('and equals the closed form <rho/(P+1)>/<P/(P+1)>', worst_c .LT. tol*scale)
  CALL check('the adiabatic electrons carry no flux surface averaged density', worst_e .LT. tol*scale)
  CALL check('<phi> is not zero: the closure was exercised', used_fsa)
  CALL check('the former <rho/P> shortcut is off by more than 1 % here', worst_k .GT. 1.0e-2_xp*scale)
  CALL check('phi is zero at the origin', ALL(phi(iky0,ikx0,:) .EQ. (0._xp,0._xp)))

  !------------------------------------------------------ ky /= 0 untouched
  CALL section('the ky /= 0 modes are the plain solve phi = rho/(P+1)')
  worst = 0._xp; scale = 0._xp
  DO ikx = 1,local_nkx
    DO iky = 1,local_nky
      IF (iky .EQ. iky0) CYCLE
      CALL polarisation(iky, ikx, P)
      CALL charge_density(iky, ikx, rho)
      DO iz = 1,local_nz
        izi = iz + ngz_o2
        worst = MAX(worst, ABS(rho(iz) - (P(iz) + 1._xp)*phi(iky,ikx,izi)), &
                           ABS(phi(iky,ikx,izi) - inv_poisson_op(iky,ikx,iz)*rho(iz)))
        scale = MAX(scale, ABS(phi(iky,ikx,izi)))
      ENDDO
    ENDDO
  ENDDO
  CALL check('the fields are not identically zero', scale .GT. 0._xp)
  CALL check('max residual of rho = (P+1) phi and of phi = inv_poisson_op rho is round off', worst .LT. tol*scale)

  CALL finish_suite

CONTAINS

  !> P(z) = sum_a q_a^2/tau_a (1 - sum_j K_j^2) from the code's kernels.
  SUBROUTINE polarisation(iky_, ikx_, P_)
    INTEGER,  INTENT(IN)  :: iky_, ikx_
    REAL(xp), INTENT(OUT) :: P_(T_NZ)
    REAL(xp) :: sumker
    INTEGER  :: iz_, izi_, ij_, iji_, ia_
    DO iz_ = 1,local_nz
      izi_ = iz_ + ngz_o2
      P_(iz_) = 0._xp
      DO ia_ = 1,local_na
        sumker = 0._xp
        DO ij_ = 1,local_nj
          iji_ = ij_ + ngj_o2
          sumker = sumker + kernel(ia_,iji_,iky_,ikx_,izi_,ieven)**2
        ENDDO
        P_(iz_) = P_(iz_) + q2_tau(ia_)*(1._xp - sumker)
      ENDDO
    ENDDO
  END SUBROUTINE polarisation

  !> rho(z) = sum_a q_a sum_j K_j N_a^{0j} from the code's kernels and moments.
  SUBROUTINE charge_density(iky_, ikx_, rho_)
    INTEGER,     INTENT(IN)  :: iky_, ikx_
    COMPLEX(xp), INTENT(OUT) :: rho_(T_NZ)
    INTEGER :: iz_, izi_, ij_, iji_, ia_
    DO iz_ = 1,local_nz
      izi_ = iz_ + ngz_o2
      rho_(iz_) = 0._xp
      DO ij_ = 1,local_nj
        iji_ = ij_ + ngj_o2
        DO ia_ = 1,local_na
          rho_(iz_) = rho_(iz_) + q(ia_)*kernel(ia_,iji_,iky_,ikx_,izi_,ieven) &
                                       *moments(ia_,ip0,iji_,iky_,ikx_,izi_,updatetlevel)
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE charge_density

  !> Random p = 0 moments on the interior j and z, zero elsewhere (ghosts, p > 0).
  SUBROUTINE random_p0_moments
    REAL(xp) :: re, im
    INTEGER  :: ia_, ij_, iky_, ikx_, iz_
    moments = 0._xp
    DO iz_ = 1,local_nz
      DO ikx_ = 1,local_nkx
        DO iky_ = 1,local_nky
          DO ij_ = 1,local_nj
            DO ia_ = 1,local_na
              CALL RANDOM_NUMBER(re); CALL RANDOM_NUMBER(im)
              moments(ia_,ip0,ij_+ngj_o2,iky_,ikx_,iz_+ngz_o2,updatetlevel) = CMPLX(re-0.5_xp, im-0.5_xp, xp)
            ENDDO
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
    WRITE(lu_in,'(A,F8.3)')  '  Lx   = ', T_LX
    WRITE(lu_in,'(A,I0)')    '  Ny   = ', T_NY
    WRITE(lu_in,'(A,F8.3)')  '  Ly   = ', T_LY
    WRITE(lu_in,'(A,I0)')    '  Nz   = ', T_NZ
    WRITE(lu_in,'(A)')       '  SG   = .false.'
    WRITE(lu_in,'(A)')       '  Nexc = 1'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&GEOMETRY'
    WRITE(lu_in,'(A)')       "  geom        = 's-alpha'"
    WRITE(lu_in,'(A)')       '  q0          = 1.4'
    WRITE(lu_in,'(A)')       '  shear       = 0.8'
    WRITE(lu_in,'(A,F8.3)')  '  eps         = ', T_EPS
    WRITE(lu_in,'(A)')       "  parallel_bc = 'dirichlet'"
    WRITE(lu_in,'(A)')       '  Npol        = 1'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&DIAGNOSTICS'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&MODEL'
    WRITE(lu_in,'(A)')       "  LINEARITY = 'linear'"
    WRITE(lu_in,'(A)')       '  Na      = 1'
    WRITE(lu_in,'(A)')       '  ADIAB_E = .true.'
    WRITE(lu_in,'(A)')       '  beta    = 0.0'
    WRITE(lu_in,'(A)')       '  nu      = 0.0'
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
    REWIND(lu_in)
  END SUBROUTINE write_input_deck

END PROGRAM test_adiab_e
