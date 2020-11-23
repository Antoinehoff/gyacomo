subroutine auxval
  !   Set auxiliary values, at beginning of simulation

  USE basic
  USE grid
  USE array
  USE fourier, ONLY: initialize_FFT
  use prec_const
  IMPLICIT NONE

  INTEGER :: irows,irowe, irow, icol
  IF (my_id .EQ. 0) WRITE(*,*) '=== Set auxiliary values ==='

  CALL set_krgrid

  CALL set_kzgrid ! Distributed dimension

  CALL set_pj

  CALL memory ! Allocate memory for global arrays

  CALL initialize_FFT ! intialization of FFTW plans

END SUBROUTINE auxval
