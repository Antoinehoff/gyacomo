SUBROUTINE ppexit
  !   Exit parallel environment

  USE basic
  USE fourier, ONLY : finalize_plans

  use prec_const
  IMPLICIT NONE


  CALL finalize_plans
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  CALL mpi_finalize(ierr)

END SUBROUTINE ppexit
