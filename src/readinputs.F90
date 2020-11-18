SUBROUTINE readinputs
  ! Additional data specific for a new run

  USE grid,             ONLY: grid_readinputs
  USE diagnostics_par,  ONLY: output_par_readinputs
  USE model,            ONLY: model_readinputs
  USE initial_par,      ONLY: initial_readinputs
  USE time_integration, ONLY: time_integration_readinputs

  USE prec_const
  IMPLICIT NONE

  !WRITE(*,'(a/)') '=== Define additional data for a new run ==='

  ! The input file must be in the same order as the following read routines commands (eg spacegrid then output then model)

  ! Load grid data from input file
  CALL grid_readinputs

  ! Load diagnostic options from input file
  CALL output_par_readinputs

  ! Load model parameters from input file
  CALL model_readinputs

  ! Load initial condition parameters from input file
  CALL initial_readinputs

  ! Load parameters for time integration from input file
  CALL time_integration_readinputs

END SUBROUTINE readinputs
