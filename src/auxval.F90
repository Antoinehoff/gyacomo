subroutine auxval
  !   Set auxiliary values, at beginning of simulation

  USE basic
  USE grid
  USE array
  ! use mumps_bsplines, only: putrow_mumps => putrow, updtmat_mumps => updtmat, factor_mumps => factor, bsolve_mumps => bsolve ! from BSPLINES
  use prec_const
  IMPLICIT NONE

  INTEGER :: irows,irowe, irow, icol
  IF (my_id .EQ. 0) WRITE(*,*) '=== Set auxiliary values ==='

  CALL set_krgrid
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)

  CALL set_kzgrid
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)

  CALL set_pj
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)

  CALL memory ! Allocate memory for global arrays
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)

END SUBROUTINE auxval
