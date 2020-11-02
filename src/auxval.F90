subroutine auxval
  !   Set auxiliary values, at beginning of simulation

  USE basic
  USE grid
  USE array
  ! use mumps_bsplines, only: putrow_mumps => putrow, updtmat_mumps => updtmat, factor_mumps => factor, bsolve_mumps => bsolve ! from BSPLINES
  use prec_const
  IMPLICIT NONE

  INTEGER :: irows,irowe, irow, icol
  WRITE(*,*) '=== Set auxiliary values ==='

  CALL set_krgrid
  CALL set_kzgrid
  CALL set_pj

  CALL memory ! Allocate memory for global arrays
END SUBROUTINE auxval
