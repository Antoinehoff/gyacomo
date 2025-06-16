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

  PUBLIC :: set_updatetlevel, time_integration_readinputs, time_integration_outputinputs, &
            adaptive_time_scheme_setup, adaptive_set_error, adaptive_must_recompute_step, &
            set_updatetlevel_rhs, adaptive_accuracy_ratio, estimate_error, dt_min

CONTAINS

  SUBROUTINE time_integration_readinputs
    ! Read the input parameters
    USE prec_const
    USE basic, ONLY : lu_in
    IMPLICIT NONE
    NAMELIST /TIME_INTEGRATION/ numerical_scheme
    NAMELIST /TIME_INTEGRATION/ adaptive_safety, adaptive_error_atol, adaptive_error_rtol, &
                               time_scheme_is_adaptive, dt_min, dt_max

    READ(lu_in,time_integration)
    CALL set_numerical_scheme
  END SUBROUTINE time_integration_readinputs
  
  SUBROUTINE time_integration_outputinputs(fid)
    USE futils, ONLY: attach, creatd
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: fid
    CHARACTER(len=256)  :: str
    
    WRITE(str,'(a)') '/data/input/time_integration'
    CALL creatd(fid, 0,(/0/),TRIM(str),'Time Integration Parameters Input')
    CALL attach(fid, TRIM(str),       "numerical_scheme", numerical_scheme)
    CALL attach(fid, TRIM(str), "time_scheme_is_adaptive", time_scheme_is_adaptive)
    CALL attach(fid, TRIM(str),        "adaptive_safety", adaptive_safety)
    CALL attach(fid, TRIM(str),     "adaptive_error_atol", adaptive_error_atol)
    CALL attach(fid, TRIM(str),     "adaptive_error_rtol", adaptive_error_rtol)
    CALL attach(fid, TRIM(str),                 "dt_min", dt_min)
    CALL attach(fid, TRIM(str),                 "dt_max", dt_max)
  END SUBROUTINE time_integration_outputinputs
  
  SUBROUTINE set_updatetlevel(new_updatetlevel)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: new_updatetlevel
    updatetlevel = new_updatetlevel
  END SUBROUTINE set_updatetlevel
  
  SUBROUTINE set_updatetlevel_rhs
    IMPLICIT NONE
    SELECT CASE (numerical_scheme)
    CASE ('RKC2')
      IF (updatetlevel == 1) THEN
        updatetlevel_rhs = 1
      ELSE
        updatetlevel_rhs = 2
      END IF
    CASE DEFAULT
      updatetlevel_rhs = updatetlevel
    END SELECT
  END SUBROUTINE set_updatetlevel_rhs
  
  !> Computes the next time step for adaptive time schemes
  !>
  !> Provided a relative stepping error, defined by `tolerance/error` where
  !> `tolerance` is the admissible error and `error` the estimated error of the
  !> current step, compute an optimal `dt` for the next step.
  !>
  !> @param[in] dt Current time step
  !> @param[in] errscale Normalized error of the current step
  !> @returns Optimal `dt` for the next step
  PURE FUNCTION adaptive_compute_next_dt(dt, errscale) RESULT(next_dt)
    REAL(xp), INTENT(IN) :: dt, errscale
    REAL(xp), PARAMETER :: MINSCALE=.1_xp, MAXSCALE=5._xp
    REAL(xp) :: next_dt, error_exponent
    
    ! The parameters in this routine seem to work "fine". Attempts to be too
    ! ambitious, especially with `adaptive_safety`, usually increase the rate
    ! at which time steps are discarded!
    IF(errscale > 1._xp) THEN
      error_exponent = 1 / REAL(adaptive_error_estimator_order + 1, xp)
    ELSE
      error_exponent = 1 / REAL(adaptive_error_estimator_order, xp)
    END IF
    
    next_dt = dt * MAX(MINSCALE, MIN(MAXSCALE, adaptive_safety * errscale**error_exponent))
    
    ! Enforce min and max dt limits
    next_dt = MAX(dt_min, MIN(dt_max, next_dt))
  END FUNCTION adaptive_compute_next_dt
  
  !> Setup for adaptive time schemes
  !>
  !> Prepares the next full embedded Runge-Kutta step (not substep) by adjusting
  !> basic::dt based on the previous time step, along with extra cleanup in case
  !> the previous step was discarded.
  SUBROUTINE adaptive_time_scheme_setup()
    USE basic, ONLY: dt, change_dt, speak
    USE parallel, ONLY: my_id
    IMPLICIT NONE
    REAL(xp) :: old_dt
    
    IF (.NOT. time_scheme_is_adaptive) RETURN
    
    old_dt = dt
    CALL change_dt(adaptive_compute_next_dt(dt, adaptive_accuracy_ratio))
    
    IF(should_discard_step) THEN
      ! Note this relates to the previous step!
      IF(my_id==0 .AND. dt/=old_dt) THEN
        CALL speak("Step discarded. dt update " // TRIM(str(old_dt)) // " => " // TRIM(str(dt)), 1)
        CALL speak("errscale was " // TRIM(str(adaptive_accuracy_ratio)), 1)
      END IF
    END IF
    
    ! Setting `adaptive_accuracy_ratio` is equivalent to saying that the step
    ! has zero error, i.e. "assume the step is perfect now, and update later
    ! when actually computing the error after performing the step"
    should_discard_step = .FALSE.
    adaptive_accuracy_ratio = HUGE(adaptive_accuracy_ratio)
  END SUBROUTINE adaptive_time_scheme_setup
  
  !> Updates the normalized adaptive time scheme error
  !>
  !> This routine is useful to combine the normalized error of the full system
  !> of equation by recording the worst (highest) error, which will eventually
  !> be passed to time_integration::adaptive_compute_next_dt.
  !>
  !> @param[in] errscale Normalized time step error
  !>
  !> @sa time_integration::adaptive_compute_next_dt
  SUBROUTINE adaptive_set_error(errscale)
    IMPLICIT NONE
    REAL(xp), INTENT(IN) :: errscale
    LOGICAL :: satisfies_tolerance
    
    IF (.NOT. time_scheme_is_adaptive) RETURN
    
    satisfies_tolerance = errscale <= 1._xp
    adaptive_accuracy_ratio = MIN(adaptive_accuracy_ratio, errscale)
    should_discard_step = should_discard_step .OR. (.NOT. satisfies_tolerance)
  END SUBROUTINE adaptive_set_error
  
  !> Returns `.true.` if the current step should be discarded.
  PURE FUNCTION adaptive_must_recompute_step() RESULT(res)
    IMPLICIT NONE
    LOGICAL :: res
    res = should_discard_step .AND. time_scheme_is_adaptive
  END FUNCTION adaptive_must_recompute_step
  
  !> Set the numerical scheme and configure appropriate parameters
  SUBROUTINE set_numerical_scheme()
    USE basic, ONLY: speak
    IMPLICIT NONE
    
    ! Basic scheme configuration
    SELECT CASE (numerical_scheme)
    CASE ('RK4')
      ntimelevel = 4
      time_scheme_is_adaptive = .FALSE.
      CALL setup_RK4()
      
    CASE ('DOPRI54')
      ntimelevel = 7  ! 5th order + 4th order for error estimation + 1 (zero-based indexing)
      time_scheme_is_adaptive = .TRUE.
      adaptive_error_estimator_order = 4
      CALL setup_DOPRI54()
      
    CASE ('RKF45')
      ntimelevel = 6  ! 4th order + 5th order for error estimation
      time_scheme_is_adaptive = .TRUE.
      adaptive_error_estimator_order = 4
      CALL setup_RKF45()
      
    CASE DEFAULT
      CALL speak("WARNING: Unrecognized numerical_scheme '" // TRIM(numerical_scheme) // "'. Using default RK4.", 0)
      numerical_scheme = 'RK4'
      ntimelevel = 4
      time_scheme_is_adaptive = .FALSE.
      CALL setup_RK4()
    END SELECT
    
    ntimelevel_rhs = ntimelevel
  END SUBROUTINE set_numerical_scheme
  
  !> Configure coefficients for classic 4th order Runge-Kutta method
  SUBROUTINE setup_RK4()
    IMPLICIT NONE
    
    IF (ALLOCATED(A_E)) DEALLOCATE(A_E)
    IF (ALLOCATED(b_E)) DEALLOCATE(b_E)
    IF (ALLOCATED(c_E)) DEALLOCATE(c_E)
    
    ALLOCATE(A_E(4,3), b_E(4,1), c_E(4))
    
    ! A_E coefficients for RK4
    A_E = 0.0_xp
    A_E(2,1) = 0.5_xp
    A_E(3,1) = 0.0_xp
    A_E(3,2) = 0.5_xp
    A_E(4,1) = 0.0_xp
    A_E(4,2) = 0.0_xp
    A_E(4,3) = 1.0_xp
    
    ! b_E coefficients for RK4
    b_E(1,1) = 1.0_xp/6.0_xp
    b_E(2,1) = 1.0_xp/3.0_xp
    b_E(3,1) = 1.0_xp/3.0_xp
    b_E(4,1) = 1.0_xp/6.0_xp
    
    ! c_E coefficients for time-dependent RHS
    c_E(1) = 0.0_xp
    c_E(2) = 0.5_xp
    c_E(3) = 0.5_xp
    c_E(4) = 1.0_xp
  END SUBROUTINE setup_RK4
  
  !> Configure coefficients for Dormand-Prince 5(4) method (DOPRI54)
  SUBROUTINE setup_DOPRI54()
    IMPLICIT NONE
    
    IF (ALLOCATED(A_E)) DEALLOCATE(A_E)
    IF (ALLOCATED(b_E)) DEALLOCATE(b_E)
    IF (ALLOCATED(c_E)) DEALLOCATE(c_E)
    
    ALLOCATE(A_E(7,6), b_E(7,2), c_E(7))
    
    ! A_E coefficients
    A_E = 0.0_xp
    A_E(2,1) = 1.0_xp/5.0_xp
    A_E(3,1) = 3.0_xp/40.0_xp
    A_E(3,2) = 9.0_xp/40.0_xp
    A_E(4,1) = 44.0_xp/45.0_xp
    A_E(4,2) = -56.0_xp/15.0_xp
    A_E(4,3) = 32.0_xp/9.0_xp
    A_E(5,1) = 19372.0_xp/6561.0_xp
    A_E(5,2) = -25360.0_xp/2187.0_xp
    A_E(5,3) = 64448.0_xp/6561.0_xp
    A_E(5,4) = -212.0_xp/729.0_xp
    A_E(6,1) = 9017.0_xp/3168.0_xp
    A_E(6,2) = -355.0_xp/33.0_xp
    A_E(6,3) = 46732.0_xp/5247.0_xp
    A_E(6,4) = 49.0_xp/176.0_xp
    A_E(6,5) = -5103.0_xp/18656.0_xp
    A_E(7,1) = 35.0_xp/384.0_xp
    A_E(7,2) = 0.0_xp
    A_E(7,3) = 500.0_xp/1113.0_xp
    A_E(7,4) = 125.0_xp/192.0_xp
    A_E(7,5) = -2187.0_xp/6784.0_xp
    A_E(7,6) = 11.0_xp/84.0_xp
    
    ! b_E coefficients (first column: 5th order, second column: 4th order)
    b_E(:,1) = 0.0_xp
    b_E(1,1) = 35.0_xp/384.0_xp
    b_E(3,1) = 500.0_xp/1113.0_xp
    b_E(4,1) = 125.0_xp/192.0_xp
    b_E(5,1) = -2187.0_xp/6784.0_xp
    b_E(6,1) = 11.0_xp/84.0_xp
    b_E(7,1) = 0.0_xp
    
    b_E(:,2) = 0.0_xp
    b_E(1,2) = 5179.0_xp/57600.0_xp
    b_E(3,2) = 7571.0_xp/16695.0_xp
    b_E(4,2) = 393.0_xp/640.0_xp
    b_E(5,2) = -92097.0_xp/339200.0_xp
    b_E(6,2) = 187.0_xp/2100.0_xp
    b_E(7,2) = 1.0_xp/40.0_xp
    
    ! c_E coefficients for time-dependent RHS
    c_E(1) = 0.0_xp
    c_E(2) = 0.2_xp
    c_E(3) = 0.3_xp
    c_E(4) = 0.8_xp
    c_E(5) = 8.0_xp/9.0_xp
    c_E(6) = 1.0_xp
    c_E(7) = 1.0_xp
  END SUBROUTINE setup_DOPRI54
  
  !> Configure coefficients for Runge-Kutta-Fehlberg 4(5) method
  SUBROUTINE setup_RKF45()
    IMPLICIT NONE
    
    IF (ALLOCATED(A_E)) DEALLOCATE(A_E)
    IF (ALLOCATED(b_E)) DEALLOCATE(b_E)
    IF (ALLOCATED(c_E)) DEALLOCATE(c_E)
    
    ALLOCATE(A_E(6,5), b_E(6,2), c_E(6))
    
    ! A_E coefficients
    A_E = 0.0_xp
    A_E(2,1) = 1.0_xp/4.0_xp
    A_E(3,1) = 3.0_xp/32.0_xp
    A_E(3,2) = 9.0_xp/32.0_xp
    A_E(4,1) = 1932.0_xp/2197.0_xp
    A_E(4,2) = -7200.0_xp/2197.0_xp
    A_E(4,3) = 7296.0_xp/2197.0_xp
    A_E(5,1) = 439.0_xp/216.0_xp
    A_E(5,2) = -8.0_xp
    A_E(5,3) = 3680.0_xp/513.0_xp
    A_E(5,4) = -845.0_xp/4104.0_xp
    A_E(6,1) = -8.0_xp/27.0_xp
    A_E(6,2) = 2.0_xp
    A_E(6,3) = -3544.0_xp/2565.0_xp
    A_E(6,4) = 1859.0_xp/4104.0_xp
    A_E(6,5) = -11.0_xp/40.0_xp
    
    ! b_E coefficients (first column: 4th order, second column: 5th order)
    b_E(:,1) = 0.0_xp
    b_E(1,1) = 25.0_xp/216.0_xp
    b_E(3,1) = 1408.0_xp/2565.0_xp
    b_E(4,1) = 2197.0_xp/4104.0_xp
    b_E(5,1) = -1.0_xp/5.0_xp
    
    b_E(:,2) = 0.0_xp
    b_E(1,2) = 16.0_xp/135.0_xp
    b_E(3,2) = 6656.0_xp/12825.0_xp
    b_E(4,2) = 28561.0_xp/56430.0_xp
    b_E(5,2) = -9.0_xp/50.0_xp
    b_E(6,2) = 2.0_xp/55.0_xp
    
    ! c_E coefficients for time-dependent RHS
    c_E(1) = 0.0_xp
    c_E(2) = 0.25_xp
    c_E(3) = 0.375_xp
    c_E(4) = 12.0_xp/13.0_xp
    c_E(5) = 1.0_xp
    c_E(6) = 0.5_xp
  END SUBROUTINE setup_RKF45
  
  !> Calculate error estimation for adaptive time stepping
  !> This estimates the local truncation error for the current step
  SUBROUTINE estimate_error(errscale)
    USE fields, ONLY: moments
    USE grid, ONLY: local_na, local_np, local_nj, local_nky, local_nkx, local_nz, &
                   ngp, ngj, ngz
    IMPLICIT NONE
    REAL(xp), INTENT(OUT) :: errscale
    COMPLEX(xp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: high_order, low_order
    REAL(xp) :: err, tol, max_err_ratio
    INTEGER :: ia, ip, ij, iky, ikx, iz
    
    ! Only perform error estimation for adaptive schemes
    IF (.NOT. time_scheme_is_adaptive) THEN
      errscale = 0.0_xp
      RETURN
    END IF
    
    ! Allocate temporary arrays for error estimation
    ALLOCATE(high_order(local_na, local_np+ngp, local_nj+ngj, local_nky, local_nkx, local_nz+ngz))
    ALLOCATE(low_order(local_na, local_np+ngp, local_nj+ngj, local_nky, local_nkx, local_nz+ngz))
    
    ! Initialize with zeros
    high_order = 0.0_xp
    low_order = 0.0_xp
    
    ! Calculate high and low order solutions for error estimation
    SELECT CASE (numerical_scheme)
    CASE ('DOPRI54')
      ! For DOPRI54, first column in b_E is 5th order, second is 4th order
      high_order = moments(:,:,:,:,:,:,1)
      low_order = moments(:,:,:,:,:,:,1)
      
    CASE ('RKF45')
      ! For RKF45, first column in b_E is 4th order, second is 5th order
      high_order = moments(:,:,:,:,:,:,1)
      low_order = moments(:,:,:,:,:,:,1)
      
    CASE DEFAULT
      ! For non-adaptive schemes, return zero error
      errscale = 0.0_xp
      DEALLOCATE(high_order, low_order)
      RETURN
    END SELECT
    
    ! Calculate the maximum relative error
    max_err_ratio = 0.0_xp
    DO iz = 1, local_nz
      DO ikx = 1, local_nkx
        DO iky = 1, local_nky
          DO ij = 1, local_nj
            DO ip = 1, local_np
              DO ia = 1, local_na
                ! Calculate absolute error
                err = ABS(high_order(ia,ip,ij,iky,ikx,iz) - low_order(ia,ip,ij,iky,ikx,iz))
                
                ! Calculate tolerance (atol + rtol * |value|)
                tol = adaptive_error_atol + adaptive_error_rtol * ABS(high_order(ia,ip,ij,iky,ikx,iz))
                
                ! Update maximum error ratio if this one is larger
                IF (tol > 0.0_xp) THEN
                  max_err_ratio = MAX(max_err_ratio, err / tol)
                END IF
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
    
    ! Compute normalized error scale (1/error used to compute next dt)
    IF (max_err_ratio > 0.0_xp) THEN
      errscale = 1.0_xp / max_err_ratio
    ELSE
      errscale = HUGE(errscale)
    END IF
    
    DEALLOCATE(high_order, low_order)
  END SUBROUTINE estimate_error
  
  FUNCTION str(k) RESULT(str_)
    !   "Convert a real to string."
    REAL(xp), INTENT(IN) :: k
    CHARACTER(LEN=32):: str_
    WRITE(str_, "(G12.4)") k
    str_ = ADJUSTL(str_)
  END FUNCTION str
  
END MODULE time_integration
