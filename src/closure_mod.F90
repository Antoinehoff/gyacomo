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
  NAMELIST /CLOSURE/ hierarchy_closure, dmax, nonlinear_closure, nmax
  READ(lu_in,closure)
END SUBROUTINE

SUBROUTINE set_closure_model
  USE basic, ONLY: speak
  USE grid,  ONLY: local_np, ngp, local_nj, ngj, parray, jarray,&
                   pmax, jmax
  USE CLA,   ONLY: set_monomial_trunc_coeff
  IMPLICIT NONE
  INTEGER :: ip,ij
  ! adapt the dmax if it is set <0
  IF(dmax .LT. 0) &
    dmax = MIN(pmax,2*jmax+1)
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
  CASE('monomial')
    DO ip = 1,local_np+ngp
      DO ij = 1, local_nj+ngj
        evolve_mom(ip,ij) = ((parray(ip).GE.0) .AND. (jarray(ij).GE.0)) &
                      .AND. (parray(ip).LE.pmax) .AND. (jarray(ij).LE.jmax)
      ENDDO
    ENDDO
  CASE('max_degree')
    IF (dmax .GT. (pmax+2*jmax)) &
      ERROR STOP "dmax is higher than the maximal moments degree available"
    DO ip = 1,local_np+ngp
      DO ij = 1, local_nj+ngj
        evolve_mom(ip,ij) = ((parray(ip).GE.0) .AND. (jarray(ij).GE.0)) &
                      .AND. (parray(ip)+2*jarray(ij) .LE. dmax)
        ! evolve_mom(ip,ij) = (parray(ip).GE.0) .AND. (jarray(ij).GE.0)
        ! evolve_mom(ip,ij) = ((parray(ip)+2*jarray(ij)) .LE. dmax)
        ! evolve_mom(ip,ij) = evolve_mom(ip,ij) .AND. ((parray(ip)+2*jarray(ij)) .LE. dmax)
      ENDDO
    ENDDO  
  CASE DEFAULT
    ERROR STOP "closure scheme not recognized (avail: truncation,max_degree,monomial)"
  END SELECT

  ! If monomial truncation, setup the coefficients required
  SELECT CASE(hierarchy_closure)
  CASE('monomial')
    CALL set_monomial_trunc_coeff(pmax,jmax)
  END SELECT

  ! Set the nonlinear closure scheme (truncation of sum over n in Sapj)
  ALLOCATE(nmaxarray(local_nj))
  SELECT CASE(nonlinear_closure)
  CASE('truncation','monomial')
    IF(nmax .LT. 0) THEN
      CALL speak("Set nonlinear truncation to anti Laguerre aliasing")
      DO ij = 1,local_nj
        nmaxarray(ij) = jmax - jarray(ij+ngj/2)
      ENDDO
    ELSE
      DO ij = 1,local_nj
        nmaxarray(ij) = MIN(nmax,jmax)
      ENDDO
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
  USE grid,       ONLY: local_nj,ngj,local_np,ngp,local_na,deltap
  USE fields,     ONLY: moments
  USE time_integration, ONLY: updatetlevel
  USE CLA, ONLY: c_Pp1, c_Pp2, c_Jp1
  USE parallel, ONLY: num_procs_p
  IMPLICIT NONE
  INTEGER ::ij,ip,ia, ipi, iji, ipmax, ijmax, iq, ik, iqi, iki
  ! Set to zero all non evolving moments
  DO ij = 1, local_nj+ngj
    DO ip = 1,local_np+ngp
      DO ia = 1,local_na
        IF(.NOT. evolve_mom(ip,ij))&
          moments(ia,ip,ij,:,:,:,updatetlevel) = 0._xp
      ENDDO
    ENDDO
  ENDDO  
  SELECT CASE (hierarchy_closure)
    CASE('truncation','max_degree')
      ! do nothing
    CASE('monomial')
      IF(num_procs_p .GT. 1) ERROR STOP "STOP: monomial closure is not parallelized in p"
      ! 
      ipmax = local_np+ngp/2
      ijmax = local_nj+ngj/2
      DO ia = 1,local_na
        IF(deltap .EQ. 1) THEN ! We solve for every Hermite
          !! P+1 truncation : NaP+1,j = sum_{q=0}^P c^{P+1}_q Nqj
          DO ij = 1, local_nj
            iji = ij + ngj/2
            moments(ia,ipmax+1,iji,:,:,:,updatetlevel) = 0._xp
            DO iq = 1, local_np
              iqi = iq + ngp/2
              moments(ia,ipmax+1,iji,:,:,:,updatetlevel) = moments(ia,ipmax+1,iji,:,:,:,updatetlevel) &
                + c_Pp1(iq) * moments(ia,iqi,iji,:,:,:,updatetlevel)
            ENDDO
          ENDDO
          !! P+2 truncation : NaP+2,j = sum_{q=0}^P c^{P+2}_q Nqj
          DO ij = 1, local_nj
            iji = ij + ngj/2
            moments(ia,ipmax+2,iji,:,:,:,updatetlevel) = 0._xp
            DO iq = 1, local_np
              iqi = iq + ngp/2
              moments(ia,ipmax+2,iji,:,:,:,updatetlevel) = moments(ia,ipmax+2,iji,:,:,:,updatetlevel) &
                + c_Pp2(iq) * moments(ia,iqi,iji,:,:,:,updatetlevel)
            ENDDO
          ENDDO
        ELSE
          !! P+2 truncation only : NaP+2,j = sum_{q=0}^P c^{P+2}_q Nqj
          DO ij = 1, local_nj
            iji = ij + ngj/2
            moments(ia,ipmax+1,iji,:,:,:,updatetlevel) = 0._xp
            DO iq = 1, local_np
              iqi = iq + ngp/2
              moments(ia,ipmax+1,iji,:,:,:,updatetlevel) = moments(ia,ipmax+1,iji,:,:,:,updatetlevel) &
                + c_Pp2(2*(iq-1)+1) * moments(ia,iqi,iji,:,:,:,updatetlevel)
            ENDDO
          ENDDO
        ENDIF
        !! J+1 truncation : Nap,J+1 = sum_{k=0}^J c^{J+1}_k Npk
        DO ip = 1, local_np
          ipi = ip + ngp/2
          moments(ia,ipi,ijmax+1,:,:,:,updatetlevel) = 0._xp
          DO ik = 1, local_nj
            iki = ik + ngj/2
            moments(ia,iki,ijmax+1,:,:,:,updatetlevel) = moments(ia,iki,ijmax+1,:,:,:,updatetlevel) &
              + c_Jp1(ik) * moments(ia,ipi,iki,:,:,:,updatetlevel)
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
