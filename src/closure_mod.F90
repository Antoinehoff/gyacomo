module closure
! Contains the routines to define closures
IMPLICIT NONE

PUBLIC :: apply_closure_model

CONTAINS

! Positive Oob indices are approximated with a model
SUBROUTINE apply_closure_model
  USE prec_const, ONLY: dp
  USE model,      ONLY: CLOS
  USE grid,       ONLY: local_nj,ngj, jarray,&
                        local_np,ngp, parray, dmax
  USE fields,     ONLY: moments
  USE time_integration, ONLY: updatetlevel
  IMPLICIT NONE
  INTEGER ::ij,ip,ia
  IF (CLOS .EQ. 0) THEN
    ! zero truncation, An+1=0 for n+1>nmax only
    CALL ghosts_upper_truncation
  ELSEIF (CLOS .EQ. 1) THEN
    ! truncation at highest fully represented kinetic moment
    ! e.g. Dmax = 3 means
    ! only Napj s.t. p+2j <= 3 are evolved
    ! -> (p,j) allowed are (0,0),(1,0),(0,1),(2,0),(1,1),(3,0)
    ! =>> Dmax = min(Pmax,2*Jmax+1)
    j: DO ij = 1,local_nj+ngj
    p: DO ip = 1,local_np+ngp
      IF ( parray(ip)+2*jarray(ij) .GT. dmax) THEN
        moments(ia,ip,ij,:,:,:,updatetlevel) = 0._dp
      ENDIF
    ENDDO p
    ENDDO j
  ELSE
    ERROR STOP '>> ERROR << Closure scheme not found '
  ENDIF
  CALL ghosts_lower_truncation
END SUBROUTINE apply_closure_model

! Positive Oob indices are approximated with a model
SUBROUTINE ghosts_upper_truncation
  USE prec_const, ONLY: dp
  USE grid,       ONLY: local_np,ngp,local_pmax, pmax,&
                        local_nj,ngj,local_jmax, jmax
  USE fields,           ONLY: moments
  USE time_integration, ONLY: updatetlevel
  IMPLICIT NONE
  INTEGER ::ig
  ! zero truncation, An+1=0 for n+1>nmax
    ! applies only for the processes that evolve the highest moment degree
    IF(local_jmax .GE. Jmax) THEN
      DO ig = 1,ngj/2
        moments(:,:,local_nj+ngj/2+ig,:,:,:,updatetlevel) = 0._dp
      ENDDO
    ENDIF
    ! applies only for the process that has largest p index
    IF(local_pmax .GE. Pmax) THEN
      DO ig = 1,ngp/2
        moments(:,local_np+ngp/2+ig,:,:,:,:,updatetlevel) = 0._dp
      ENDDO
    ENDIF
END SUBROUTINE ghosts_upper_truncation

! Negative OoB indices are 0
SUBROUTINE ghosts_lower_truncation
  USE prec_const, ONLY: dp
  USE grid,       ONLY: ngp,ngj,local_pmin,local_jmin
  USE fields,           ONLY: moments
  USE time_integration, ONLY: updatetlevel
  IMPLICIT NONE
  INTEGER :: ig
! zero truncation, An=0 for n<0
    IF(local_jmin .EQ. 0) THEN
      DO ig  = 1,ngj/2
        moments(:,:,ig,:,:,:,updatetlevel) = 0._dp
      ENDDO
    ENDIF
    ! applies only for the process that has lowest p index
    IF(local_pmin .EQ. 0) THEN
      DO ig  = 1,ngp/2
        moments(:,ig,:,:,:,:,updatetlevel) = 0._dp
      ENDDO
    ENDIF

END SUBROUTINE ghosts_lower_truncation

END module closure
