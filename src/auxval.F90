subroutine auxval
  !   Set auxiliary values, at beginning of simulation

  USE basic
  USE fourier_grid
  USE array
  ! use mumps_bsplines, only: putrow_mumps => putrow, updtmat_mumps => updtmat, factor_mumps => factor, bsolve_mumps => bsolve ! from BSPLINES
  use prec_const
  IMPLICIT NONE
  
  INTEGER :: irows,irowe, irow, icol
  WRITE(*,*) '=== Set auxiliary values ==='

  call set_krgrid
  call set_kzgrid
  call set_pj

  CALL memory ! Allocate memory for global arrays
END SUBROUTINE auxval
