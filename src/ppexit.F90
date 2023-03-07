SUBROUTINE ppexit
  !   Exit parallel environment

  USE basic
  USE fourier, ONLY : finalize_plans

  use prec_const
  IMPLICIT NONE
  INTEGER :: ierr

  CALL finalize_plans
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  CALL MPI_FINALIZE(ierr)

END SUBROUTINE ppexit
