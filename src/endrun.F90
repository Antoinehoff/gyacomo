SUBROUTINE endrun
  !   Terminate the run

  USE basic
  use prec_const
  IMPLICIT NONE


  IF( nlend ) THEN
     !----------------------------------------------------------------------
     !              1.   Normal end of run
     WRITE(*,'(/a)') '   Normal exit'

     !----------------------------------------------------------------------
     !              2.   Abnormal exit
  ELSE
     WRITE(*,'(/a)') '   Abnormal exit'
  END IF


END SUBROUTINE endrun
