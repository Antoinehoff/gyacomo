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
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)

  CALL set_krgrid
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)

  CALL set_kzgrid ! Distributed dimension
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)

  CALL set_pj
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)

  CALL memory ! Allocate memory for global arrays
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)

END SUBROUTINE auxval
