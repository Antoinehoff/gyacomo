SUBROUTINE tesend
  !   Test for run completion

  USE basic

  use prec_const
  IMPLICIT NONE

  !________________________________________________________________________________
  !                   1.  Some processors had set nlend
  IF( nlend ) THEN
    WRITE(*,'(/a)') 'rhs are NaN/Inf'
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
  CALL cpu_time(finish)
  nlend = 1.1*(finish-start) .GT. maxruntime
  IF ( nlend ) THEN
     IF (my_id .EQ. 0) WRITE(*,'(/a)') 'Max run time reached'
     RETURN
  END IF
  !
  !
END SUBROUTINE tesend
