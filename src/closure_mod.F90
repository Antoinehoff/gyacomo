module closure
! Contains the routines to define closures
USE basic
USE model,  ONLY: CLOS, tau_e, tau_i, q_e, q_i, nu, KIN_E
USE grid
USE array,  ONLY: kernel_e,  kernel_i
USE fields, ONLY: moments_e, moments_i
USE time_integration, ONLY: updatetlevel
IMPLICIT NONE

PUBLIC :: apply_closure_model, ghosts_truncation

CONTAINS

! Positive Oob indices are approximated with a model
SUBROUTINE apply_closure_model
  IMPLICIT NONE

  CALL cpu_time(t0_clos)
  IF (CLOS .EQ. 0) THEN
    ! zero truncation, An+1=0 for n+1>nmax only
    CALL ghosts_truncation


  ELSEIF (CLOS .EQ. 1) THEN
    ! Truncation at highest fully represented kinetic moment
    ! e.g. Dmax = 3 means
    ! all Napj s.t. p+2j <= 3
    ! -> (p,j) allowed are (0,0),(1,0),(0,1),(2,0),(1,1),(3,0)
    ! =>> Dmax is Pmax, condition is p+2j<=Pmax
    DO ikx = ikxs,ikxe
      DO iky = ikys,ikye
        DO iz = izs,ize
          IF(KIN_E) THEN
          DO ip = ipsg_e,ipeg_e
            DO ij = ijsg_e,ijeg_e
              IF ( parray_e(ip)+2*jarray_e(ip) .GT. dmaxe) &
              moments_e(ip,ij,ikx,iky,iz,updatetlevel) = 0._dp
            ENDDO
          ENDDO
          ENDIF
          DO ip = ipsg_i,ipeg_i
            DO ij = ijsg_i,ijeg_i
              IF ( parray_i(ip)+2*jarray_i(ip) .GT. dmaxi) &
              moments_i(ip,ij,ikx,iky,iz,updatetlevel) = 0._dp
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    ! + ghosts truncation
    CALL ghosts_truncation
  ELSE
    if(my_id .EQ. 0) write(*,*) '! Closure scheme not found !'

  ENDIF

  CALL cpu_time(t1_clos)
  tc_clos = tc_clos + (t1_clos - t0_clos)
END SUBROUTINE apply_closure_model

! Positive Oob indices are approximated with a model
SUBROUTINE ghosts_truncation
  IMPLICIT NONE

! zero truncation, An+1=0 for n+1>nmax
    DO ikx = ikxs,ikxe
      DO iky = ikys,ikye
        DO iz = izs,ize
          IF(KIN_E) THEN
          DO ip = ipsg_e,ipeg_e
            moments_e(ip,ijsg_e,ikx,iky,iz,updatetlevel) = 0._dp
            moments_e(ip,ijeg_e,ikx,iky,iz,updatetlevel) = 0._dp
          ENDDO
          DO ij = ijsg_e,ijeg_e
            moments_e(ipsg_e  ,ij,ikx,iky,iz,updatetlevel) = 0._dp
            moments_e(ipeg_e  ,ij,ikx,iky,iz,updatetlevel) = 0._dp
            IF(deltape .EQ. 1) THEN ! Must truncate the second stencil
            moments_e(ipsg_e+1,ij,ikx,iky,iz,updatetlevel) = 0._dp
            moments_e(ipeg_e-1,ij,ikx,iky,iz,updatetlevel) = 0._dp
            ENDIF
          ENDDO
          kernel_e(ijsg_e,ikx,iky,iz)      = 0._dp
          kernel_e(ijeg_e,ikx,iky,iz)      = 0._dp
          ENDIF
          DO ip = ipsg_i,ipeg_i
            moments_i(ip,ijsg_i,ikx,iky,iz,updatetlevel) = 0._dp
            moments_i(ip,ijeg_i,ikx,iky,iz,updatetlevel) = 0._dp
          ENDDO
          DO ij = ijsg_i,ijeg_i
            moments_i(ipsg_i  ,ij,ikx,iky,iz,updatetlevel) = 0._dp
            moments_i(ipeg_i  ,ij,ikx,iky,iz,updatetlevel) = 0._dp
            IF(deltapi .EQ. 1) THEN ! Must truncate the second stencil
            moments_i(ipsg_i+1,ij,ikx,iky,iz,updatetlevel) = 0._dp
            moments_i(ipeg_i-1,ij,ikx,iky,iz,updatetlevel) = 0._dp
            ENDIF
          ENDDO
          kernel_i(ijsg_i,ikx,iky,iz)      = 0._dp
          kernel_i(ijeg_i,ikx,iky,iz)      = 0._dp
        ENDDO
      ENDDO
    ENDDO

END SUBROUTINE ghosts_truncation

END module closure
