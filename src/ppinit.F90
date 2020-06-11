SUBROUTINE ppinit
  !   Parallel environment

  USE basic

  use prec_const
  IMPLICIT NONE
  
  INTEGER :: version_prov=-1

  ! CALL mpi_init(ierr)
  ! CALL MPI_INIT_THREAD(MPI_THREAD_SINGLE,version_prov,ierr)
  CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,version_prov,ierr)

END SUBROUTINE ppinit
