MODULE time_integration

  USE prec_const

  IMPLICIT NONE
  PRIVATE

  INTEGER, PUBLIC, PROTECTED :: ntimelevel=4 ! Total number of time levels required by the numerical scheme
  INTEGER, PUBLIC, PROTECTED :: updatetlevel ! Current time level to be updated

  real(xp),PUBLIC,PROTECTED,DIMENSION(:,:),ALLOCATABLE :: A_E
  real(xp),PUBLIC,PROTECTED,DIMENSION(:,:),ALLOCATABLE :: b_E
  real(xp),PUBLIC,PROTECTED,DIMENSION(:),ALLOCATABLE :: c_E !Coeff for Expl/Implic time integration in case of time dependent RHS (i.e. dy/dt = f(y,t)) see Baptiste Frei CSE Rapport 06/17

  character(len=10),PUBLIC,PROTECTED :: numerical_scheme='RK4'

  !!---- start: This part is taken from GBS adaptive time stepping by L. Stenger
  ! For more information check gbs.epfl.ch
  logical, public, protected :: time_scheme_is_adaptive = .false.
  logical, public, protected :: should_discard_step = .false.
  real(xp) :: adaptive_accuracy_ratio = 1._xp
  integer  :: adaptive_error_estimator_order = -1
  real(xp) :: dt_min = 1e-8_xp          !< Minimum allowable timestep for adaptive time schemes
  real(xp) :: dt_max = 1e1_xp           !< Maximum allowable timestep for adaptive time schemes
  real(xp) :: adaptive_safety = 0.9_xp  !< Safety factor of the adaptive timestep controller
  !> @{
  !> @brief Error tolerance for adaptive time schemes
  !> @details
  !> Typically, an adaptive time scheme will provide an error estimation for a
  !> given time step "witdh". The allowable error will be defined as
  !> `atol + rtol * abs(f)` (where `f` the quantity which is being integrated).
  real(xp), public, protected :: adaptive_error_atol=1e-6_xp, adaptive_error_rtol=1e-3_xp
  !> @}
  integer, public, protected :: ntimelevel_rhs   = 4 !< number of time levels required for the rhs
  integer, public, protected :: updatetlevel_rhs = 1 !< time level to be updated for rhs
  !!---- end

  PUBLIC :: set_updatetlevel, time_integration_readinputs, time_integration_outputinputs

CONTAINS

  SUBROUTINE time_integration_readinputs
    ! Read the input parameters
    USE prec_const
    USE basic, ONLY : lu_in
    IMPLICIT NONE
    NAMELIST /TIME_INTEGRATION/ numerical_scheme
    namelist /TIME_INTEGRATION/ adaptive_safety, adaptive_error_atol, adaptive_error_rtol

    READ(lu_in,time_integration)
    CALL set_numerical_scheme
  END SUBROUTINE time_integration_readinputs


  SUBROUTINE time_integration_outputinputs(fid)
    ! Write the input parameters to the results_xx.h5 file
    USE futils, ONLY: attach, creatd
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fid
    CHARACTER(len=256)  :: str
    WRITE(str,'(a)') '/data/input/time_integration'
    CALL creatd(fid, 0,(/0/),TRIM(str),'Time Integration Input')
    CALL attach(fid, TRIM(str), "numerical_scheme", numerical_scheme)
  END SUBROUTINE time_integration_outputinputs

  SUBROUTINE set_numerical_scheme
    ! Initialize Butcher coefficient of set_numerical_scheme
    use parallel, ONLY: my_id
    IMPLICIT NONE
    SELECT CASE (numerical_scheme)
    ! Order 2 methods
    CASE ('RK2')
      CALL RK2
    ! Order 3 methods
    CASE ('RK3')
      CALL RK3
    CASE ('SSP_RK3')
      CALL SSP_RK3
    ! Adaptative scheme
    CASE ('ODE23')
      time_scheme_is_adaptive = .true.
      CALL ode23
    ! Order 4 methods
    CASE ('RK4')
      CALL RK4
    ! Order 5 methods
    CASE ('DOPRI5')
      time_scheme_is_adaptive = .true.
      CALL DOPRI5
    CASE DEFAULT
       IF (my_id .EQ. 0) WRITE(*,*) 'Cannot initialize time integration scheme. Name invalid.'
    END SELECT
    IF (my_id .EQ. 0) WRITE(*,*) " Time integration with ", numerical_scheme
  END SUBROUTINE set_numerical_scheme

  !!! second order time schemes
  SUBROUTINE RK2
    ! Butcher coeff for clasical RK2 (Heun's)
    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 2
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep,1,1)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 2
    c_E(1)   = 0.0_xp
    c_E(2)   = 1.0_xp
    b_E(1,1) = 1._xp/2._xp
    b_E(2,1) = 1._xp/2._xp
    A_E(2,1) = 1._xp
  END SUBROUTINE RK2

  !!! third order time schemes
  SUBROUTINE RK3
    ! Butcher coeff for classical RK3
    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 3
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep,1,1)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 3
    c_E(1)   = 0.0_xp
    c_E(2)   = 1.0_xp/2.0_xp
    c_E(3)   = 1.0_xp
    b_E(1,1) = 1._xp/6._xp
    b_E(2,1) = 2._xp/3._xp
    b_E(3,1) = 1._xp/6._xp
    A_E(2,1) = 1.0_xp/2.0_xp
    A_E(3,1) = -1._xp
    A_E(3,2) = 2._xp
  END SUBROUTINE RK3

  SUBROUTINE SSP_RK3
    ! Butcher coeff for strong stability  preserving RK3
    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 3
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep,1,1)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 3
    c_E(1)   = 0.0_xp
    c_E(2)   = 1.0_xp
    c_E(3)   = 1.0_xp/2.0_xp
    b_E(1,1) = 1._xp/6._xp
    b_E(2,1) = 1._xp/6._xp
    b_E(3,1) = 2._xp/3._xp
    A_E(2,1) = 1._xp
    A_E(3,1) = 1._xp/4._xp
    A_E(3,2) = 1._xp/4._xp
  END SUBROUTINE SSP_RK3

  !!! fourth order time schemes
  SUBROUTINE RK4
    ! Butcher coeff for RK4 (default)
    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 4
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep,1,1)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 4
    c_E(1)   = 0.0_xp
    c_E(2)   = 1.0_xp/2.0_xp
    c_E(3)   = 1.0_xp/2.0_xp
    c_E(4)   = 1.0_xp
    b_E(1,1) = 1.0_xp/6.0_xp
    b_E(2,1) = 1.0_xp/3.0_xp
    b_E(3,1) = 1.0_xp/3.0_xp
    b_E(4,1) = 1.0_xp/6.0_xp
    A_E(2,1) = 1.0_xp/2.0_xp
    A_E(3,2) = 1.0_xp/2.0_xp
    A_E(4,3) = 1.0_xp
  END SUBROUTINE RK4

  !!! fifth order time schemes
  SUBROUTINE DOPRI5
    ! Butcher coeff for DOPRI5 --> Stiffness detection
    ! DOPRI5 used for stiffness detection.
    ! 5 order method/7 stages
    USE basic
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep =7
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep,1,2)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)

    adaptive_error_estimator_order = 4
    ntimelevel = 7
    ntimelevel_rhs = 7

    c_E = 0._xp
    c_E(2) = 1.0_xp/5.0_xp
    c_E(3) = 3.0_xp /10.0_xp
    c_E(4) = 4.0_xp/5.0_xp
    c_E(5) = 8.0_xp/9.0_xp
    c_E(6) = 1.0_xp
    c_E(7) = 1.0_xp

    A_E = 0._xp
    A_E(2,1) = 1.0_xp/5.0_xp
    A_E(3,1) = 3.0_xp/40.0_xp
    A_E(3,2) = 9.0_xp/40.0_xp
    A_E(4,1) = 44.0_xp/45.0_xp
    A_E(4,2) = -56.0_xp/15.0_xp
    A_E(4,3) = 32.0_xp/9.0_xp
    A_E(5,1 ) = 19372.0_xp/6561.0_xp
    A_E(5,2) = -25360.0_xp/2187.0_xp
    A_E(5,3) = 64448.0_xp/6561.0_xp
    A_E(5,4) = -212.0_xp/729.0_xp
    A_E(6,1) = 9017.0_xp/3168.0_xp
    A_E(6,2)= -355.0_xp/33.0_xp
    A_E(6,3) = 46732.0_xp/5247.0_xp
    A_E(6,4) = 49.0_xp/176.0_xp
    A_E(6,5) = -5103.0_xp/18656.0_xp
    A_E(7,1) = 35.0_xp/384.0_xp
    A_E(7,3) = 500.0_xp/1113.0_xp
    A_E(7,4) = 125.0_xp/192.0_xp
    A_E(7,5) = -2187.0_xp/6784.0_xp
    A_E(7,6) = 11.0_xp/84.0_xp

    b_E = 0._xp
    b_E(1, 1) = 35._xp/384._xp
    b_E(2, 1) = 0._xp
    b_E(3, 1) = 500._xp/1113._xp
    b_E(4, 1) = 125._xp/192._xp
    b_E(5, 1) = -2187._xp/6784._xp
    b_E(6, 1) = 11._xp/84._xp
    b_E(7, 1) = 0._xp
    b_E(1, 2) = 5179._xp/57600._xp
    b_E(2, 2) = 0._xp
    b_E(3, 2) = 7571._xp/16695._xp
    b_E(4, 2) = 393._xp/640._xp
    b_E(5, 2) = -92097._xp/339200._xp
    b_E(6, 2) = 187._xp/2100._xp
    b_E(7, 2) = 1._xp/40._xp
  END SUBROUTINE DOPRI5
  !!-------------------------------------------------------------------------
  !!---- This part is taken from GBS adaptive time stepping by L. Stenger ---
  !> Updates the time step
  !>
  !> If the new time step in not within the closed interval `[dt_min, dt_max]`,
  !> the simulation will be terminated
  !>
  !> @param[in] new_dt Desired time step value
  subroutine set_dt(new_dt)
    use basic, only: dt, nlend, change_dt
    use parallel, ONLY: comm_p, my_id
    real(xp), intent(in) :: new_dt
    integer :: ierr
    logical :: mlend
    mlend = .false.
    if(dt<=dt_min .or. dt>=dt_max) then
      mlend = .true.
      if(my_id == 0) then
        print "(A, 3(G12.5,A))", "The adaptive time integration scheme tried to set the time step to ", dt, &
          ", which is outside of the allowed range of [", dt_min, ", ", dt_max, &
          "]." // new_line("A") // "The simulation will be terminated."
      end if
    else
      CALL change_dt(new_dt)
    end if
    call mpi_allreduce(mlend, nlend, 1, MPI_LOGICAL, MPI_LOR, comm_p, ierr)
  end subroutine set_dt
  !> Computes the next time step for adaptive time schemes
  !>
  !> Provided a relative stepping error, defined by `tolerance/error` where
  !> `tolerance` is the admissible error and `error` the estimated error of the
  !> current step, compute an optimal `dt` for the next step.
  !>
  !> @param[in] dt Current time step
  !> @param[in] errscale Normalized error of the current step
  !> @returns Optimal `dt` for the next step
  pure function adaptive_compute_next_dt(dt, errscale) result(next_dt)
    real(xp), intent(in) :: dt, errscale
    real(xp), parameter :: MINSCALE=.1_xp, MAXSCALE=5._xp
    real(xp) :: next_dt, error_exponent
    ! The parameters in this routine seem to work "fine". Attempts to be too
    ! ambitious, especially with `adaptive_safety`, usually increase the rate
    ! at which time steps are discarded!
    if(errscale > 1._xp) then
        error_exponent = 1 / real(adaptive_error_estimator_order + 1, xp)
    else
        error_exponent = 1 / real(adaptive_error_estimator_order, xp)
    end if
    next_dt =  dt * max(MINSCALE, min(MAXSCALE, adaptive_safety * errscale**error_exponent))
  end function adaptive_compute_next_dt   

  !> Setup for adaptive time schemes
  !>
  !> Prepares the next full embedded Runge-Kutta step (not substep) by adjusting
  !> basic::dt based on the previous time step, along with extra cleanup in case
  !> the previous step was discarded.
  subroutine adaptive_time_scheme_setup()
    use basic, only: dt, str, speak
    real(xp) :: old_dt

    old_dt = dt
    call set_dt(adaptive_compute_next_dt(dt, adaptive_accuracy_ratio))

    if(should_discard_step) then
        ! Note this relates to the previous step!
        if(dt/=old_dt) CALL speak("step discarded. dt update "//str(old_dt)//" => "//str(dt))
        if(dt/=old_dt) CALL speak("errscale was "//str(adaptive_accuracy_ratio))
    end if

    ! Setting `adaptive_accuracy_ratio` is equivalent to saying that the step
    ! has zero error, i.e. "assume the step is perfect now, and update later
    ! when actually computing the error after performing the step"
    should_discard_step = .false.
    adaptive_accuracy_ratio = huge(adaptive_accuracy_ratio)
  end subroutine adaptive_time_scheme_setup

  !> Updates the normalized adaptive time scheme error
  !>
  !> This routine is useful to combine the normalized error of the full system
  !> of equation by recording the worst (highest) error, which will eventually
  !> be passed to time_integration::adaptive_compute_next_dt.
  !>
  !> @param[in] errscale Normalized time step error
  !>
  !> @sa time_integration::adaptive_compute_next_dt
  subroutine adaptive_set_error(errscale)
    real(xp), intent(in) :: errscale
    logical :: satisfies_tolerance
    satisfies_tolerance = errscale > 1._xp
    adaptive_accuracy_ratio = min(adaptive_accuracy_ratio, errscale)
    should_discard_step = should_discard_step .or. (.not. satisfies_tolerance)
  end subroutine adaptive_set_error

  !> Returns `.true.` if the current step should be discarded.
  pure function adaptive_must_recompute_step() result(res)
  logical :: res
  res = should_discard_step
  end function adaptive_must_recompute_step

  subroutine set_updatetlevel(new_updatetlevel)
  integer, intent(in) :: new_updatetlevel
  updatetlevel = new_updatetlevel
  end subroutine set_updatetlevel

  subroutine set_updatetlevel_rhs
  select case (numerical_scheme)
  case ('rkc2')
  if (updatetlevel == 1)then
    updatetlevel_rhs = 1
    else
    updatetlevel_rhs = 2
    end if
  case default
    updatetlevel_rhs = updatetlevel
  end select
  end subroutine set_updatetlevel_rhs

  !> Sets up Butcher coefficients for the Bogacki-Shampine method (orders 3 and 2)
  subroutine ode23
    use basic
    integer,parameter :: nbstep = 4
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep, 1, 2)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)

    adaptive_error_estimator_order = 2
    ntimelevel = 4
    ntimelevel_rhs = 4

    c_E = 0._xp
    c_E(2) = 1.0_xp/2.0_xp
    c_E(3) = 3.0_xp /4.0_xp
    c_E(4) = 1._xp

    A_E = 0._xp
    A_E(2,1) = 1.0_xp/2.0_xp
    A_E(3,2) = 3.0_xp/4.0_xp
    A_E(4,1) = 2.0_xp/9.0_xp
    A_E(4,2) = 1.0_xp/3.0_xp
    A_E(4,3) = 4.0_xp/9.0_xp

    b_E = 0._xp
    b_E(1, 1) = 2._xp/9._xp
    b_E(2, 1) = 1._xp/3._xp
    b_E(3, 1) = 4._xp/9._xp
    b_E(4, 1) = 0._xp
    b_E(1, 2) = 7._xp/24._xp
    b_E(2, 2) = 1._xp/4._xp
    b_E(3, 2) = 1._xp/3._xp
    b_E(4, 2) = 1._xp/8._xp
  end subroutine ode23
  !!-------------------------------------------------------------------------
  !!-------------------------------------------------------------------------

END MODULE time_integration
