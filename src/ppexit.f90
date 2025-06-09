SUBROUTINE ppexit
  !   Exit parallel environment

  USE basic
  USE fourier, ONLY : finalize_plans
  USE model,   ONLY : LINEARITY

  use prec_const
  IMPLICIT NONE
  INTEGER :: ierr

  IF (LINEARITY .EQ. 'nonlinear') &
   CALL finalize_plans
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  CALL MPI_FINALIZE(ierr)

END SUBROUTINE ppexit
