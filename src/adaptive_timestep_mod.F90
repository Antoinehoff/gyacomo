!! This module is a compilation of the routines implemented in GBS by L. Stenger
!! for adaptive time stepping.
MODULE adaptive_timestep
    USE prec_const, ONLY: xp

    IMPLICIT NONE

    logical, public, protected :: time_scheme_is_adaptive = .false.
    logical, public, protected :: should_discard_step = .false.
    real(xp) :: adaptive_accuracy_ratio = 1._dp
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
    real(xp), public, protected :: adaptive_error_atol=1e-6_dp, adaptive_error_rtol=1e-3_dp
    !> @}
    integer, public, protected :: ntimelevel_rhs=4   !< number of time levels required for the rhs
  

    CONTAINS


        SUBROUTINE adaptive_timestep_readinputs
            ! Read the input parameters
            USE prec_const
            USE basic, ONLY : lu_in
            IMPLICIT NONE
            namelist /TIME_INTEGRATION/ adaptive_safety, adaptive_error_atol, adaptive_error_rtol

            READ(lu_in,time_integration)
            CALL set_numerical_scheme
        END SUBROUTINE adaptive_timestep_readinputs

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
            if(errscale > 1._dp) then
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
            use basic, only: cstepP, dt, meW, stepP
            real(xp) :: old_dt

            old_dt = dt
            call set_dt(adaptive_compute_next_dt(dt, adaptive_accuracy_ratio))

            if(should_discard_step) then
                ! Note this relates to the previous step!
                if(meW==0 .and. dt/=old_dt) print *, "step discarded. dt update ", old_dt, " => ", dt
                if(meW==0 .and. dt/=old_dt) print *, "errscale was ", adaptive_accuracy_ratio
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

END  MODULE adaptive_timestep