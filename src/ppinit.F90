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
  
END SUBROUTINE ppinit
