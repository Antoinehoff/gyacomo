SUBROUTINE tesend
  !   Test for run completion

  USE basic

  use prec_const
  IMPLICIT NONE

  !________________________________________________________________________________
  !                   1.  Some processors had set nlend
  IF( nlend ) THEN
    WRITE(*,'(/a)') 'rhs are NaN/Inf'  
    WRITE(*,*) 'Run terminated at cstep=',cstep
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
     WRITE(*,'(/a)') 'TMAX reached'
     RETURN
  END IF
  ! 
  !
END SUBROUTINE tesend
