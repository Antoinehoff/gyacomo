SUBROUTINE ppinit
  !   Parallel environment

  USE basic
  USE model, ONLY : NON_LIN
  use prec_const
  IMPLICIT NONE

  INTEGER :: version_prov=-1

  CALL MPI_INIT(ierr)
  ! CALL MPI_INIT_THREAD(MPI_THREAD_SINGLE,version_prov,ierr)
  ! CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,version_prov,ierr)

  CALL MPI_COMM_RANK (MPI_COMM_WORLD,     my_id, ierr)
  CALL MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
  CALL MPI_COMM_GROUP(MPI_COMM_WORLD,grp_world,ierr)
  ! CALL MPI_COMM_CREATE_GROUP(MPI_COMM_WORLD,grp_world,0,comm_self,ierr)

  ! IF (NON_LIN .AND. (num_procs .GT. 1)) THEN
  !   CALL fftw_mpi_init
  ! ENDIF

END SUBROUTINE ppinit
