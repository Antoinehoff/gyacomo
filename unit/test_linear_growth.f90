PROGRAM test_linear_growth
  ! Pins down the on-the-fly eigenvalue measurement of linear_growth_mod.
  !
  ! The module answers the question a linear run is launched to answer, so a
  ! silent error in it is expensive: a wrong gamma looks exactly like a right
  ! one. There is no reference to compare against at run time either, which is
  ! why the estimator is driven here with a field whose eigenvalue is known by
  ! construction,
  !
  !   phi(t) = A exp((g - i w) t),
  !
  ! and the measurement is required to return that g and that w.
  !
  ! The second thing checked here is subtler and is the reason the test exists
  ! at all. The module renormalises the state to stop a growing mode from
  ! overflowing, and the measurement is a ratio of the field at two times. If a
  ! rescale lands between those two times and only one side of the ratio is
  ! scaled, the measured growth rate picks up a spurious log(factor)/dt which is
  ! not obviously wrong when you look at it: it is a smooth, plausible number.
  ! So the invariance of gamma and omega across a rescale is asserted directly.
  !
  ! Single rank, so every local size equals its total.
  USE prec_const, ONLY: xp, imagu, PI
  USE basic,      ONLY: lu_in, dt, time, increase_time, change_dt
  USE parallel,   ONLY: ppinit
  USE grid,       ONLY: total_nky, local_nky, ngz_o2, kyarray
  USE fields,     ONLY: phi
  USE model,      ONLY: LINEARITY
  USE linear_growth
  USE test_harness
  IMPLICIT NONE

  ! Ny = 8 gives 5 ky rows, so ky = 0 sits next to four rows that must all be
  ! measured: enough to catch a spectrum that is only right in its first entry.
  INTEGER,  PARAMETER :: T_PMAX = 2, T_JMAX = 1
  INTEGER,  PARAMETER :: T_NX = 4, T_NY = 8, T_NZ = 4
  REAL(xp), PARAMETER :: T_LX = 60._xp, T_LY = 60._xp
  INTEGER,  PARAMETER :: T_NA = 1
  ! Growth rates and frequencies of the driven modes, one pair per ky row.
  ! Mixed signs so that a dropped minus sign cannot pass, and the phase turn per
  ! measurement is kept under pi so the complex logarithm stays single valued.
  REAL(xp), PARAMETER :: G_WANT(5) = (/ 0._xp,  0.8_xp, -0.3_xp,  1.7_xp, 0.05_xp /)
  REAL(xp), PARAMETER :: W_WANT(5) = (/ 0._xp, -1.2_xp,  0.7_xp, -0.4_xp, 2.00_xp /)
  ! Interval between two measurements. The estimator is exact for a pure
  ! exponential whatever this is, so the only constraint is the branch cut of
  ! the complex logarithm: the phase must turn by less than pi between two
  ! measurements, i.e. max|w|*dt_m = 2.0*0.5 = 1.0 < pi. readinputs does not set
  ! dt (only basic_data does, and this test does not call it), so the step is
  ! imposed here rather than read.
  INTEGER,  PARAMETER :: N_SUB   = 5
  REAL(xp), PARAMETER :: T_DT    = 0.1_xp
  ! Change of a drifting rate between two measurements: under the deck's
  ! lin_tol = 1e-8 relative to rates of order one, while lin_nconv = 3 of them
  ! add up to 1.2e-8, over it. Well clear of round off either way.
  REAL(xp), PARAMETER :: DRIFT_STEP = 0.4e-8_xp

  REAL(xp) :: tol, dt_m, amp_before, amp_after, phase
  COMPLEX(xp) :: phi_unscaled
  INTEGER  :: iky, im
  LOGICAL  :: ok

  CALL ppinit
  CALL start_suite('linear_growth')
  tol = tol_of(1.0e5_xp)

  CALL write_input_deck
  CALL readinputs
  CALL auxval
  CLOSE(lu_in)
  CALL change_dt(T_DT)

  dt_m = N_SUB*dt

  !------------------------------------------------------- configuration
  CALL section('the run was configured the way the test assumes')
  CALL check('LINEARITY is linear',   TRIM(LINEARITY) .EQ. 'linear')
  CALL check('LIN_GROWTH was read',   LIN_GROWTH)
  CALL check_int('5 ky rows',         total_nky, T_NY/2+1)
  CALL check_int('single rank in ky', local_nky, total_nky)
  CALL check_real('ky(1) is the zonal row', kyarray(1), 0._xp, tol)
  ! Without this the complex logarithm wraps and the frequencies come back
  ! aliased, which is a property of the driving used here and not of the code.
  phase = MAXVAL(ABS(W_WANT))*dt_m
  CALL check('the driving turns by less than pi per measurement', phase .LT. PI)

  !------------------------------------------------------- the estimator
  ! Two measurements are needed before anything is reported: the first only
  ! records the reference field.
  CALL set_phi
  CALL lin_measure
  CALL advance(N_SUB)
  CALL set_phi
  CALL lin_measure
  CALL advance(N_SUB)
  CALL set_phi
  CALL lin_measure

  CALL section('a known exp((g - i w)t) is recovered per ky')
  DO iky = 2,total_nky
    CALL check_real('gamma at ky row '//istr(iky), gamma_ky(iky), G_WANT(iky), tol)
    CALL check_real('omega at ky row '//istr(iky), omega_ky(iky), W_WANT(iky), tol)
  ENDDO

  ! The zonal row is skipped rather than measured: it is not an eigenmode here
  ! and its projection would be 0/0.
  CALL section('the zonal row is left alone')
  CALL check_real('gamma stays zero at ky = 0', gamma_ky(1), 0._xp, tol)
  CALL check_real('omega stays zero at ky = 0', omega_ky(1), 0._xp, tol)

  ! Only two values have been produced so far, fewer than the lin_nconv = 3 the
  ! deck asks for, so the run must still want to continue however steady they are.
  CALL section('convergence does not latch before lin_nconv measurements')
  CALL check('has not latched yet', .NOT. lin_converged)

  !------------------------------------------------------- rescaling
  ! Let the fastest mode grow until it leaves the band on its own, the way it
  ! would in a run, rather than planting a large number by hand. Measuring on
  ! every pass is also what a run does, and is what keeps the phase turn per
  ! measurement inside the branch the complex logarithm is single valued on.
  CALL section('rescaling brings the amplitude back into the band')
  DO im = 1,200
    CALL advance(N_SUB)
    CALL set_phi
    IF (MAXVAL(ABS(phi)) .GT. lin_ceiling) EXIT
    CALL lin_measure
  ENDDO
  amp_before   = MAXVAL(ABS(phi))
  ! The amplitude this mode would have carried had nothing ever been divided out.
  phi_unscaled = phi(2,1,1+ngz_o2)*EXP(lin_lnscale)
  CALL check('the mode grew past the ceiling', amp_before .GT. lin_ceiling)
  CALL lin_rescale
  amp_after    = MAXVAL(ABS(phi))
  CALL check('amplitude is back under the ceiling', amp_after .LE. lin_ceiling)
  CALL check('amplitude is still above the floor',  amp_after .GE. lin_floor)
  ! What is written to disk times EXP(lin_lnscale) is what the run would have
  ! held without any renormalisation; that is the whole contract of lin_lnscale.
  CALL check_cmpx('lin_lnscale undoes the rescaling', &
                  phi(2,1,1+ngz_o2)*EXP(lin_lnscale), phi_unscaled, tol)

  CALL section('the eigenvalue is blind to the rescaling')
  ! phi_prev was rescaled by the same factor inside lin_rescale, so this
  ! measurement spans the rescale event and must still see the true rate.
  CALL advance(N_SUB)
  CALL set_phi
  CALL lin_measure
  DO iky = 2,total_nky
    CALL check_real('gamma across a rescale, ky row '//istr(iky), &
                    gamma_ky(iky), G_WANT(iky), tol)
    CALL check_real('omega across a rescale, ky row '//istr(iky), &
                    omega_ky(iky), W_WANT(iky), tol)
  ENDDO

  !------------------------------------------------------- slow drifts
  ! A steady drift whose step between two measurements is under lin_tol, but
  ! which adds up to more than lin_tol across the lin_nconv+1 values of the
  ! window, is a rate that is still moving: the run must not stop on it. This
  ! is checked on gamma and on omega separately, since both have to settle.
  CALL section('a slow drift in gamma alone does not latch')
  DO im = 1,3*(lin_nconv+1)
    CALL advance(N_SUB)
    CALL set_phi(ramp_g=DRIFT_STEP/dt_m)
    CALL lin_rescale
    CALL lin_measure
  ENDDO
  CALL check('drifting gamma keeps the run going', .NOT. lin_converged)

  CALL section('a slow drift in omega alone does not latch')
  DO im = 1,3*(lin_nconv+1)
    CALL advance(N_SUB)
    CALL set_phi(ramp_w=DRIFT_STEP/dt_m)
    CALL lin_rescale
    CALL lin_measure
  ENDDO
  CALL check('drifting omega keeps the run going', .NOT. lin_converged)

  !------------------------------------------------------- convergence
  ! A field whose eigenvalue never moves must end up being declared settled,
  ! including across the renormalisations that happen along the way.
  CALL section('convergence latches once the rate holds still')
  ! Twice the window: the first measurement after the ramps are switched off
  ! sees the quadratic term vanish at once, a one-off jump that has to leave
  ! the window before the values can agree again.
  DO im = 1,2*(lin_nconv+1)
    CALL advance(N_SUB)
    CALL set_phi
    CALL lin_rescale       ! as the main loop does, every step
    CALL lin_measure
  ENDDO
  CALL check('latched once the rate held still', lin_converged)

  CALL finish_suite

CONTAINS

  !> Impose the field a converged linear run would hold at the current time,
  !! phi = EXP((g - i w) t) per ky row, with the zonal row left at zero.
  !!
  !! The EXP(-lin_lnscale) factor is what makes this an honest stand-in for the
  !! solver. Once the state has been renormalised, the code does not go back to
  !! the true amplitude: it keeps integrating the scaled field, and lin_lnscale
  !! is exactly the log of what has been divided out so far. Dropping the factor
  !! here would re-inflate the field behind the module's back and show up as a
  !! spurious jump in the growth rate, which is an artefact of the test rather
  !! than a fault in the code.
  !!
  !! ramp_g and ramp_w add a constant d(gamma)/dt and d(omega)/dt, through a
  !! quadratic term in the exponent, so that the measured rate walks steadily.
  !! They are relative rates, scaled by each row's own gamma or omega, so that
  !! every row sits just under lin_tol per measurement: a row far under or far
  !! over it would decide the outcome on its own and hide what is being tested.
  SUBROUTINE set_phi(ramp_g, ramp_w)
    REAL(xp), INTENT(IN), OPTIONAL :: ramp_g, ramp_w
    REAL(xp) :: ag, aw
    INTEGER  :: jky
    ag = 0._xp
    aw = 0._xp
    IF (PRESENT(ramp_g)) ag = ramp_g
    IF (PRESENT(ramp_w)) aw = ramp_w
    phi = (0._xp,0._xp)
    DO jky = 2,local_nky
      phi(jky,:,:) = EXP((G_WANT(jky) - imagu*W_WANT(jky))*time &
                         + 0.5_xp*(ag*ABS(G_WANT(jky)) - imagu*aw*ABS(W_WANT(jky)))*time**2 &
                         - lin_lnscale)
    ENDDO
  END SUBROUTINE set_phi

  !> Step the clock without touching the fields: the test supplies those itself.
  SUBROUTINE advance(nstep)
    INTEGER, INTENT(IN) :: nstep
    INTEGER :: is
    DO is = 1,nstep
      CALL increase_time
    ENDDO
  END SUBROUTINE advance

  !> Small integer to string, for assertion labels.
  CHARACTER(len=8) FUNCTION istr(i_)
    INTEGER, INTENT(IN) :: i_
    WRITE(istr,'(I0)') i_
  END FUNCTION istr

  !> Feed the namelists through a scratch file so nothing is left on disk. Only
  !! the entries this test depends on are set; the rest keep their defaults.
  SUBROUTINE write_input_deck
    OPEN(UNIT=lu_in, STATUS='SCRATCH', ACTION='READWRITE', FORM='FORMATTED')
    WRITE(lu_in,'(A)')      '&BASIC'
    WRITE(lu_in,'(A)')      '  dt   = 0.1'
    WRITE(lu_in,'(A)')      '  tmax = 1000'
    WRITE(lu_in,'(A)')      '/'
    WRITE(lu_in,'(A)')      '&GRID'
    WRITE(lu_in,'(A,I0)')   '  pmax = ', T_PMAX
    WRITE(lu_in,'(A,I0)')   '  jmax = ', T_JMAX
    WRITE(lu_in,'(A,I0)')   '  Nx   = ', T_NX
    WRITE(lu_in,'(A,F8.2)') '  Lx   = ', T_LX
    WRITE(lu_in,'(A,I0)')   '  Ny   = ', T_NY
    WRITE(lu_in,'(A,F8.2)') '  Ly   = ', T_LY
    WRITE(lu_in,'(A,I0)')   '  Nz   = ', T_NZ
    WRITE(lu_in,'(A)')      '  SG   = .false.'
    WRITE(lu_in,'(A)')      '  Nexc = 1'
    WRITE(lu_in,'(A)')      '/'
    ! z-pinch keeps the geometry out of the way; the estimator does not use it.
    WRITE(lu_in,'(A)')      '&GEOMETRY'
    WRITE(lu_in,'(A)')      "  geom = 'z-pinch'"
    WRITE(lu_in,'(A)')      '/'
    WRITE(lu_in,'(A)')      '&DIAGNOSTICS'
    WRITE(lu_in,'(A)')      '/'
    WRITE(lu_in,'(A)')      '&MODEL'
    WRITE(lu_in,'(A)')      "  LINEARITY = 'linear'"
    WRITE(lu_in,'(A,I0)')   '  Na      = ', T_NA
    WRITE(lu_in,'(A)')      '  ADIAB_E = .true.'
    WRITE(lu_in,'(A)')      '  beta    = 0.0'
    WRITE(lu_in,'(A)')      '  nu      = 0.0'
    WRITE(lu_in,'(A)')      '/'
    WRITE(lu_in,'(A)')      '&CLOSURE'
    WRITE(lu_in,'(A)')      "  hierarchy_closure = 'truncation'"
    WRITE(lu_in,'(A)')      "  nonlinear_closure = 'truncation'"
    WRITE(lu_in,'(A)')      '  nmax = -1'
    WRITE(lu_in,'(A)')      '/'
    WRITE(lu_in,'(A)')      '&SPECIES'
    WRITE(lu_in,'(A)')      "  name_  = 'ions'"
    WRITE(lu_in,'(A)')      '  tau_   = 1.0'
    WRITE(lu_in,'(A)')      '  sigma_ = 1.0'
    WRITE(lu_in,'(A)')      '  q_     = 1.0'
    WRITE(lu_in,'(A)')      '  k_N_   = 2.0'
    WRITE(lu_in,'(A)')      '  k_T_   = 0.4'
    WRITE(lu_in,'(A)')      '/'
    WRITE(lu_in,'(A)')      '&COLLISION'
    WRITE(lu_in,'(A)')      "  collision_model = 'DG'"
    WRITE(lu_in,'(A)')      '  GK_CO           = .true.'
    WRITE(lu_in,'(A)')      '/'
    WRITE(lu_in,'(A)')      '&INITIAL'
    WRITE(lu_in,'(A)')      "  INIT_OPT = 'phi'"
    WRITE(lu_in,'(A)')      '/'
    WRITE(lu_in,'(A)')      '&TIME_INTEGRATION'
    WRITE(lu_in,'(A)')      "  numerical_scheme = 'RK4'"
    WRITE(lu_in,'(A)')      '/'
    ! The group under test. A tight tolerance is fine because the driven field
    ! is a pure exponential: successive measurements agree to round off.
    WRITE(lu_in,'(A)')      '&LINEAR_GROWTH'
    WRITE(lu_in,'(A)')      '  LIN_GROWTH  = .true.'
    WRITE(lu_in,'(A)')      '  lin_tol     = 1.0e-8'
    WRITE(lu_in,'(A)')      '  lin_nconv   = 3'
    WRITE(lu_in,'(A)')      '  lin_ceiling = 1.0e3'
    WRITE(lu_in,'(A)')      '  lin_floor   = 1.0e-3'
    WRITE(lu_in,'(A)')      '/'
    REWIND(lu_in)
  END SUBROUTINE write_input_deck

END PROGRAM test_linear_growth
