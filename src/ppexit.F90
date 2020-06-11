SUBROUTINE ppexit
  !   Exit parallel environment

  USE basic

  use prec_const
  IMPLICIT NONE
  

  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  CALL mpi_finalize(ierr)
  
END SUBROUTINE ppexit
