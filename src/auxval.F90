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
  CALL set_kpgrid

  CALL memory ! Allocate memory for global arrays

  !! Display parallel settings
  DO i_ = 0,num_procs-1
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    IF (my_id .EQ. i_) THEN
      IF (my_id .EQ. 0) WRITE(*,*) ''
      IF (my_id .EQ. 0) WRITE(*,*) '--------- Parallel environement ----------'
      IF (my_id .EQ. 0) WRITE(*,'(A9,I3,A10,I3,A10,I3)') 'n_procs= ', num_procs, ', num_procs_p   = ', num_procs_p, ', num_procs_kr   = ', num_procs_kr
      IF (my_id .EQ. 0) WRITE(*,*) ''
      WRITE(*,'(A9,I3,A10,I3,A10,I3)')&
       'my_id  = ', my_id, ', rank_p  = ', rank_p, ', rank_kr  = ', rank_kr
       WRITE(*,'(A22,I3,A10,I3)')&
       '              ips_e = ', ips_e, ', ikrs  = ', ikrs
       WRITE(*,'(A22,I3,A10,I3)')&
       '              ipe_e = ', ipe_e, ', ikre  = ', ikre
       WRITE(*,'(A22,I3,A10,I3)')&
       '              ips_i = ', ips_i, ', ikzs  = ', ikzs
       WRITE(*,'(A22,I3,A10,I3)')&
       '              ipe_i = ', ipe_i, ', ikze  = ', ikze

      !  WRITE(*,'(A9,I3,A10,I3,A10,I3,A10,I3)')&
      !   ' ips_e = ',ips_e,', ipe_e = ',ipe_e,', ips_i = ',ips_i,', ipe_i = ',ipe_i
      ! WRITE (*,'(A9,I3,A10,I3,A10,I3,A10,I3,A10,I3)') &
      ! ' ikrs  = ', ikrs, ', ikre  = ', ikre, ', ikzs  = ', ikzs, ', ikze = ', ikze
      IF (my_id .NE. num_procs-1) WRITE (*,*) ''
      IF (my_id .EQ. num_procs-1) WRITE(*,*) '------------------------------------------'
    ENDIF
  ENDDO
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)

END SUBROUTINE auxval
