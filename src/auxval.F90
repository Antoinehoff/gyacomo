subroutine auxval
  !   Set auxiliary values, at beginning of simulation

  USE basic
  USE grid
  USE array
  USE fourier, ONLY: init_grid_distr_and_plans
  use prec_const
  IMPLICIT NONE

  INTEGER :: irows,irowe, irow, icol
  IF (my_id .EQ. 0) WRITE(*,*) '=== Set auxiliary values ==='

  CALL init_grid_distr_and_plans(Nr,Nz)

  CALL set_krgrid

  CALL set_kzgrid ! Distributed dimension

  CALL set_pj

  CALL memory ! Allocate memory for global arrays

END SUBROUTINE auxval
