subroutine auxval
  !   Set auxiliary values, at beginning of simulation

  USE basic
  USE grid
  USE array
  USE model
  USE fourier, ONLY: init_grid_distr_and_plans, alloc_local_1, alloc_local_2
  use prec_const
  IMPLICIT NONE

  INTEGER :: irows,irowe, irow, icol, i_
  IF (my_id .EQ. 0) WRITE(*,*) '=== Set auxiliary values ==='

  IF (NON_LIN) THEN
    CALL init_grid_distr_and_plans(Nr,Nz)
  ELSE
    CALL init_1Dgrid_distr
  ENDIF

  CALL set_pgrid
  CALL set_jgrid

  CALL set_krgrid ! MPI Distributed dimension
  CALL set_kzgrid

  CALL memory ! Allocate memory for global arrays

  DO i_ = 0,num_procs-1
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    IF (my_id .EQ. i_) THEN
      WRITE (*,'(I2,A9,I3,A8,I3,A8,I3,A8,I3,A15,I6)') &
        i_,': ikrs = ', ikrs, ' ikre = ', ikre, 'ikzs = ', ikzs, ' ikze = ', ikze
      WRITE (*,'(A14,I4,A10,I4,A15,I6)') &
        '  local_nkr = ',local_nkr,' offset = ',local_nkr_offset,' alloc_local = ', alloc_local_1+alloc_local_2
    ENDIF
  ENDDO

END SUBROUTINE auxval
