subroutine auxval
  !   Set auxiliary values, at beginning of simulation

  USE basic
  USE grid
  USE array
  USE fourier, ONLY: init_grid_distr_and_plans, alloc_local_1, alloc_local_2
  use prec_const
  IMPLICIT NONE

  INTEGER :: irows,irowe, irow, icol, i_
  IF (my_id .EQ. 0) WRITE(*,*) '=== Set auxiliary values ==='

  CALL init_grid_distr_and_plans(Nr,Nz)

  CALL set_krgrid ! MPI Distributed dimension

  CALL set_kzgrid

  CALL set_pj

  CALL memory ! Allocate memory for global arrays

  DO i_ = 0,num_procs-1
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    IF (my_id .EQ. i_) WRITE (*,'(I2,A9,I3,A8,I3,A8,I3,A8,I3,A15,I6)') &
    i_,': ikrs = ', ikrs, ' ikre = ', ikre, 'ikzs = ', ikzs, ' ikze = ', ikze, &
    'alloc_local = ', alloc_local_1+alloc_local_2
  ENDDO

END SUBROUTINE auxval
