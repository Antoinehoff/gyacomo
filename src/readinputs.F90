SUBROUTINE readinputs
  ! Additional data specific for a new run

  USE grid,             ONLY: grid_readinputs
  USE diagnostics,      ONLY: diag_readinputs
  USE collision,        ONLY: collision_readinputs
  USE model,            ONLY: model_readinputs
  USE species,          ONLY: species_readinputs
  USE initial,          ONLY: initial_readinputs
  USE time_integration, ONLY: time_integration_readinputs
  USE geometry,         ONLY: geometry_readinputs
  USE closure,          ONLY: closure_readinputs

  USE prec_const
  IMPLICIT NONE

  !WRITE(*,'(a/)') '=== Define additional data for a new run ==='

  ! The input file must be in the same order as the following read routines commands (eg spacegrid then output then model)

  ! Load grid data from input file
  CALL grid_readinputs

  ! Load geometry parameters
  CALL geometry_readinputs

  ! Load diagnostic options from input file
  CALL diag_readinputs

  ! Load model parameters from input file
  CALL model_readinputs

  ! Load parameters for moment closure scheme
  CALL closure_readinputs
  
  ! Load model parameters from input file
  CALL species_readinputs

  ! Load collision parameters from input file
  CALL collision_readinputs

  ! Load initial condition parameters from input file
  CALL initial_readinputs

  ! Load parameters for time integration from input file
  CALL time_integration_readinputs

END SUBROUTINE readinputs
