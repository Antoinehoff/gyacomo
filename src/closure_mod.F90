module closure
! Contains the routines to define closures
IMPLICIT NONE
! Input
CHARACTER(len=32),  PUBLIC, PROTECTED :: hierarchy_closure = 'truncation'! closure for the moment hierarchy
INTEGER,            PUBLIC, PROTECTED :: dmax              = -1          ! max evolved degree moment
CHARACTER(len=32),  PUBLIC, PROTECTED :: nonlinear_closure = 'truncation'! nonlinear truncation method
INTEGER,            PUBLIC, PROTECTED :: nmax              = 0           ! upperbound of the nonlinear sum over n
!  Attributes
LOGICAL,DIMENSION(:,:), ALLOCATABLE, PUBLIC, PROTECTED :: evolve_mom     ! array that sets if a moment has to be evolved or not
INTEGER,DIMENSION(:),   ALLOCATABLE, PUBLIC, PROTECTED :: nmaxarray      ! upperbound of the nonlinear sum over n (depend on j)
 
PUBLIC :: closure_readinputs, set_closure_model, apply_closure_model

CONTAINS

SUBROUTINE closure_readinputs
  USE basic, ONLY: lu_in
  IMPLICIT NONE
  NAMELIST /CLOSURE_PAR/ hierarchy_closure, dmax, nonlinear_closure, nmax
  READ(lu_in,closure_par)
END SUBROUTINE

SUBROUTINE set_closure_model
  USE basic, ONLY: speak
  USE grid, ONLY: local_np, ngp, local_nj, ngj, parray, jarray,&
                  pmax, jmax
  IMPLICIT NONE
  INTEGER :: ip,ij
  ! adapt the dmax if it is set <0
  IF(dmax .LT. 0) THEN
    dmax = MIN(pmax,2*jmax+1)
  ELSEIF(dmax .GT. (pmax+2*jmax)) THEN
    ERROR STOP "dmax is higher than the maximal moments degree available"
  ENDIF
  ! set the evolve mom array
  ALLOCATE(evolve_mom(local_np+ngp,local_nj+ngj))
  SELECT CASE(hierarchy_closure)
  CASE('truncation')
    DO ip = 1,local_np+ngp
      DO ij = 1, local_nj+ngj
        evolve_mom(ip,ij) = ((parray(ip).GE.0) .AND. (jarray(ij).GE.0)) &
                      .AND. (parray(ip).LE.pmax) .AND. (jarray(ij).LE.jmax)
      ENDDO
    ENDDO
  CASE('max_degree')
    DO ip = 1,local_np+ngp
      DO ij = 1, local_nj+ngj
          evolve_mom(ip,ij) = ((parray(ip).GE.0) .AND. (jarray(ij).GE.0)) &
                        .AND. (parray(ip)+2*jarray(ij) .LE. dmax)
      ENDDO
    ENDDO  
  CASE DEFAULT
    ERROR STOP "closure scheme not recognized (avail: truncation,max_degree)"
  END SELECT
  ! Set the nonlinear closure scheme (truncation of sum over n in Sapj)
  ALLOCATE(nmaxarray(local_nj))
  SELECT CASE(nonlinear_closure)
  CASE('truncation')
    IF(nmax .LT. 0) THEN
      CALL speak("Set nonlinear truncation to anti Laguerre aliasing")
      DO ij = 1,local_nj
        nmaxarray(ij) = jmax - jarray(ij+ngj/2)
      ENDDO
    ELSE
      nmaxarray(:) = nmax
    ENDIF
  CASE('anti_laguerre_aliasing')
    DO ij = 1,local_nj
      nmaxarray(ij) = jmax - jarray(ij+ngj/2)
    ENDDO
  CASE('full_sum')
    nmaxarray(:) = jmax
  CASE DEFAULT
    ERROR STOP "nonlinear closure scheme not recognized (avail: truncation,anti_laguerre_aliasing,full_sum)"
  END SELECT

END SUBROUTINE set_closure_model

! Positive Oob indices are approximated with a model
SUBROUTINE apply_closure_model
  USE prec_const, ONLY: xp
  USE grid,       ONLY: local_nj,ngj,local_np,ngp,local_na
  USE fields,     ONLY: moments
  USE time_integration, ONLY: updatetlevel
  IMPLICIT NONE
  INTEGER ::ij,ip,ia
  SELECT CASE (hierarchy_closure)
    CASE('truncation','max_degree')
      DO ij = 1, local_nj+ngj
        DO ip = 1,local_np+ngp
          DO ia = 1,local_na
            IF(.NOT. evolve_mom(ip,ij))&
              moments(ia,ip,ij,:,:,:,updatetlevel) = 0._xp
          ENDDO
        ENDDO
      ENDDO  
    CASE DEFAULT
      ERROR STOP "closure scheme not recognized"
  END SELECT
END SUBROUTINE apply_closure_model


SUBROUTINE closure_outputinputs(fid)
  ! Write the input parameters to the results_xx.h5 file
  USE futils, ONLY: attach, creatd
  IMPLICIT NONE
  INTEGER, INTENT(in) :: fid
  CHARACTER(len=256)  :: str
  WRITE(str,'(a)') '/data/input/closure'
  CALL creatd(fid, 0,(/0/),TRIM(str),'Closure Input')
  CALL attach(fid, TRIM(str),"hierarchy_closure",hierarchy_closure)
  CALL attach(fid, TRIM(str),             "dmax",dmax)
  CALL attach(fid, TRIM(str),"nonlinear_closure",nonlinear_closure)
  CALL attach(fid, TRIM(str),             "nmax",nmax)
END SUBROUTINE closure_outputinputs

END module closure
