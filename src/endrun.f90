SUBROUTINE endrun
  !   Terminate the run

  USE basic
  use prec_const
  IMPLICIT NONE


  IF( nlend ) THEN
     !----------------------------------------------------------------------
     !              1.   Normal end of run
     CALL speak('Normal exit',0)

     !----------------------------------------------------------------------
     !              2.   Abnormal exit
  ELSE
     CALL speak('Abnormal exit',0)
  END IF


END SUBROUTINE endrun
