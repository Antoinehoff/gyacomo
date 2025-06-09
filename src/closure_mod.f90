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
                   pmax, jmax, ngj_o2
  USE CLA,   ONLY: set_monomial_trunc_coeff
  IMPLICIT NONE
  INTEGER :: ip,ij
  ! adapt the dmax if it is set <0
  IF(dmax .LT. 0) &
    dmax = MIN(pmax,2*jmax+1)
  ! set the evolve mom array
  ALLOCATE(evolve_mom(local_np+ngp,local_nj+ngj))
  evolve_mom = 0
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
      ENDDO
    ENDDO  
  CASE('hot_electrons')
    !! Evolve moments only up to the maximal degree MIN(pmax,2*jmax+1)
    dmax = MIN(pmax,2*jmax+1)
    DO ip = 1,local_np+ngp
      DO ij = 1, local_nj+ngj
        evolve_mom(ip,ij) = ((parray(ip).GE.0) .AND. (jarray(ij).GE.0)) &
                    .AND. (parray(ip)+2*jarray(ij) .LE. dmax)
      ENDDO
    ENDDO  
    !! Change the linear matrix coefficients to remove coupling
    ! in the higher order GM equations
    CALL asymptotic_tau_closure
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
      CALL speak("Set nonlinear truncation to anti Laguerre aliasing",2)
      DO ij = 1,local_nj
        nmaxarray(ij) = jmax - jarray(ij+ngj_o2)
      ENDDO
    ELSE
      DO ij = 1,local_nj
        nmaxarray(ij) = MIN(nmax,jmax)
      ENDDO
    ENDIF
  CASE('anti_laguerre_aliasing')
    DO ij = 1,local_nj
      nmaxarray(ij) = jmax - jarray(ij+ngj_o2)
    ENDDO
  CASE('full_sum')
    nmaxarray(:) = jmax
  CASE DEFAULT
    ERROR STOP "nonlinear closure scheme not recognized (avail: truncation,anti_laguerre_aliasing,full_sum)"
  END SELECT
END SUBROUTINE set_closure_model

! Closure based on an ordering of tau. We remove all O(tau) terms in the last moment equation
SUBROUTINE asymptotic_tau_closure
  USE array, ONLY:  xnapj, &
                    ynapp1j, ynapm1j, ynapp1jm1, ynapm1jm1,&
                    zNapm1j, zNapm1jp1, zNapm1jm1,&
                    xnapj, xnapjp1, xnapjm1,&
                    xnapp1j, xnapm1j, xnapp2j, xnapm2j
  USE grid,  ONLY:  parray, jarray, &
                    local_na, local_np, local_nj, ngj_o2, ngp_o2
  USE prec_const, ONLY: xp
  IMPLICIT NONE
  INTEGER     :: ia,ip,ij,p_int, j_int ! polynom. dagrees
  REAL(xp)    :: p_xp, j_xp
  !! these large loops are not efficient but it is done only once
  !! and it's easier to read.
  DO ia = 1,local_na
    DO ip = 1,local_np
      p_int= parray(ip+ngp_o2)   ! Hermite degree
      p_xp = REAL(p_int,xp) ! REAL of Hermite degree
      DO ij = 1,local_nj
        j_int= jarray(ij+ngj_o2)   ! Laguerre degree
        j_xp = REAL(j_int,xp) ! REAL of Laguerre degree
        !! We se to zero all effects related to tau 
        !! in the two highest degree GMs
        IF((parray(ip)+2*jarray(ij) .GE. (dmax-1))) THEN
            ! these cancel out the magn. curvature and perp. gradient coeff.
            xnapj(ia,ip,ij)= 0._xp
            xnapp2j(ia,ip) = 0._xp
            xnapm2j(ia,ip) = 0._xp
            xnapjp1(ia,ij) = 0._xp
            xnapjm1(ia,ij) = 0._xp
        ENDIF
        !! We se to zero all effects related to sqrt(tau) 
        !! in the highest degree GMs
        IF((parray(ip)+2*jarray(ij) .GE. (dmax))) THEN
          ! Mirror force terms
          ynapp1j(ia,ip,ij)   = 0._xp
          ynapm1j(ia,ip,ij)   = 0._xp
          ynapp1jm1(ia,ip,ij) = 0._xp
          ynapm1jm1(ia,ip,ij) = 0._xp
          ! Trapping terms
          zNapm1j(ia,ip,ij)   = 0._xp
          zNapm1jp1(ia,ip,ij) = 0._xp
          zNapm1jm1(ia,ip,ij) = 0._xp
          ! Landau damping
          xnapp1j(ia,ip)      = 0._xp
          xnapm1j(ia,ip)      = 0._xp
        ENDIF
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE asymptotic_tau_closure

! Positive Oob indices are approximated with a model
SUBROUTINE apply_closure_model
  USE prec_const, ONLY: xp
  USE grid,       ONLY: local_nj,ngj,local_np,ngp,local_na,deltap,ngj_o2,ngp_o2
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
    CASE('truncation','max_degree','hot_electrons')
      ! do nothing
    CASE('monomial')
      IF(num_procs_p .GT. 1) ERROR STOP "STOP: monomial closure is not parallelized in p"
      ! 
      ipmax = local_np+ngp_o2
      ijmax = local_nj+ngj_o2
      DO ia = 1,local_na
        IF(deltap .EQ. 1) THEN ! We solve for every Hermite
          !! P+1 truncation : NaP+1,j = sum_{q=0}^P c^{P+1}_q Nqj
          DO ij = 1, local_nj
            iji = ij + ngj_o2
            moments(ia,ipmax+1,iji,:,:,:,updatetlevel) = 0._xp
            DO iq = 1, local_np
              iqi = iq + ngp_o2
              moments(ia,ipmax+1,iji,:,:,:,updatetlevel) = moments(ia,ipmax+1,iji,:,:,:,updatetlevel) &
                + c_Pp1(iq) * moments(ia,iqi,iji,:,:,:,updatetlevel)
            ENDDO
          ENDDO
          !! P+2 truncation : NaP+2,j = sum_{q=0}^P c^{P+2}_q Nqj
          DO ij = 1, local_nj
            iji = ij + ngj_o2
            moments(ia,ipmax+2,iji,:,:,:,updatetlevel) = 0._xp
            DO iq = 1, local_np
              iqi = iq + ngp_o2
              moments(ia,ipmax+2,iji,:,:,:,updatetlevel) = moments(ia,ipmax+2,iji,:,:,:,updatetlevel) &
                + c_Pp2(iq) * moments(ia,iqi,iji,:,:,:,updatetlevel)
            ENDDO
          ENDDO
        ELSE
          !! P+2 truncation only : NaP+2,j = sum_{q=0}^P c^{P+2}_q Nqj
          DO ij = 1, local_nj
            iji = ij + ngj_o2
            moments(ia,ipmax+1,iji,:,:,:,updatetlevel) = 0._xp
            DO iq = 1, local_np
              iqi = iq + ngp_o2
              moments(ia,ipmax+1,iji,:,:,:,updatetlevel) = moments(ia,ipmax+1,iji,:,:,:,updatetlevel) &
                + c_Pp2(2*(iq-1)+1) * moments(ia,iqi,iji,:,:,:,updatetlevel)
            ENDDO
          ENDDO
        ENDIF
        !! J+1 truncation : Nap,J+1 = sum_{k=0}^J c^{J+1}_k Npk
        DO ip = 1, local_np
          ipi = ip + ngp_o2
          moments(ia,ipi,ijmax+1,:,:,:,updatetlevel) = 0._xp
          DO ik = 1, local_nj
            iki = ik + ngj_o2
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
