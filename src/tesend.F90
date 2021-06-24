SUBROUTINE tesend
  !   Test for run completion

  USE basic

  use prec_const
  IMPLICIT NONE
  LOGICAL :: mlend
  real    :: tnow

  !________________________________________________________________________________
  !                   1.  Some processors had set nlend
  CALL mpi_allreduce(nlend, mlend, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &
     &             ierr)
  IF( mlend ) THEN
    nlend   = .TRUE.
    crashed = .TRUE.
    IF (my_id .EQ. 0) WRITE(*,'(/a)') 'rhs are NaN/Inf'
    IF (my_id .EQ. 0) WRITE(*,*) 'Run terminated at cstep=',cstep
    RETURN
  END IF

  !________________________________________________________________________________
  !                   2.  Test on NRUN
  nlend = step .GT. nrun
  IF ( nlend ) THEN
     WRITE(*,'(/a)') 'NRUN steps done'
     RETURN
  END IF


  !________________________________________________________________________________
  !                   3.  Test on TMAX
  nlend = time .GT. tmax
  IF ( nlend ) THEN
     IF (my_id .EQ. 0) WRITE(*,'(/a)') 'TMAX reached'
     RETURN
  END IF
  !
  !

  !________________________________________________________________________________
  !                   4.  Test on rune time
  CALL cpu_time(tnow)
  mlend = (1.2*(tnow-start)) .GT. maxruntime
  CALL mpi_allreduce(mlend, nlend, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &
       &    ierr)
  IF ( nlend ) THEN
     IF(my_id.EQ.0) WRITE(*,'(/a)') 'Max run time reached'
     RETURN
  END IF
  !
  RETURN
  !
END SUBROUTINE tesend
