SUBROUTINE ppexit
  !   Exit parallel environment

  USE basic
  USE fourier, ONLY : finalize_FFT

  use prec_const
  IMPLICIT NONE


  CALL finalize_FFT
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  CALL mpi_finalize(ierr)

END SUBROUTINE ppexit
