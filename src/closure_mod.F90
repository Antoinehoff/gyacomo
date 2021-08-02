module closure
! Contains the routines to define closures
USE basic
USE model,  ONLY: CLOS, tau_e, tau_i, q_e, q_i, eta_B, nu
USE grid
USE array,  ONLY: kernel_e,  kernel_i
USE fields, ONLY: moments_e, moments_i
USE time_integration, ONLY: updatetlevel
IMPLICIT NONE

PUBLIC :: apply_closure_model

CONTAINS

! Positive Oob indices are approximated with a model
SUBROUTINE apply_closure_model
  IMPLICIT NONE

  CALL cpu_time(t0_clos)
! zero truncation, An+1=0 for n+1>nmax
  IF (CLOS .EQ. 0) THEN
    DO ikx = ikxs,ikxe
      DO iky = ikys,ikye
        DO iz = izs,ize
          DO ip = ipsg_e,ipeg_e
            moments_e(ip,ijsg_e,ikx,iky,iz,updatetlevel) = 0._dp
            moments_e(ip,ijeg_e,ikx,iky,iz,updatetlevel) = 0._dp
          ENDDO
          DO ij = ijsg_e,ijeg_e
            moments_e(ipsg_e+1,ij,ikx,iky,iz,updatetlevel) = 0._dp
            moments_e(ipsg_e  ,ij,ikx,iky,iz,updatetlevel) = 0._dp
            moments_e(ipeg_e-1,ij,ikx,iky,iz,updatetlevel) = 0._dp
            moments_e(ipeg_e  ,ij,ikx,iky,iz,updatetlevel) = 0._dp
          ENDDO
          kernel_e(ijsg_e,ikx,iky)      = 0._dp
          kernel_e(ijeg_e,ikx,iky)      = 0._dp

          DO ip = ipsg_i,ipeg_i
            moments_i(ip,ijsg_i,ikx,iky,iz,updatetlevel) = 0._dp
            moments_i(ip,ijeg_i,ikx,iky,iz,updatetlevel) = 0._dp
          ENDDO
          DO ij = ijsg_i,ijeg_i
            moments_i(ipsg_i+1,ij,ikx,iky,iz,updatetlevel) = 0._dp
            moments_i(ipsg_i  ,ij,ikx,iky,iz,updatetlevel) = 0._dp
            moments_i(ipeg_i-1,ij,ikx,iky,iz,updatetlevel) = 0._dp
            moments_i(ipeg_i  ,ij,ikx,iky,iz,updatetlevel) = 0._dp
          ENDDO
          kernel_i(ijsg_i,ikx,iky)      = 0._dp
          kernel_i(ijeg_i,ikx,iky)      = 0._dp
        ENDDO
      ENDDO
    ENDDO

  ! zero truncation, An+1=0 for n+1>nmax
  ELSEIF (CLOS .EQ. 1) THEN
    DO ikx = ikxs,ikxe
      DO iky = ikys,ikye
        DO iz = izs,ize
          DO ip = ipsg_e,ipeg_e
            moments_e(ip,ijsg_e,ikx,iky,iz,updatetlevel) = 0._dp
            moments_e(ip,ijeg_e,ikx,iky,iz,updatetlevel) = 0._dp
          ENDDO
          DO ij = ijsg_e,ijeg_e
            moments_e(ipsg_e+1,ij,ikx,iky,iz,updatetlevel) = 0._dp
            moments_e(ipsg_e  ,ij,ikx,iky,iz,updatetlevel) = 0._dp
            moments_e(ipeg_e-1,ij,ikx,iky,iz,updatetlevel) = 0._dp
            moments_e(ipeg_e  ,ij,ikx,iky,iz,updatetlevel) = 0._dp
          ENDDO
          kernel_e(ijsg_e,ikx,iky)      = 0._dp
          kernel_e(ijeg_e,ikx,iky)      = 0._dp

          DO ip = ipsg_i,ipeg_i
            moments_i(ip,ijsg_i,ikx,iky,iz,updatetlevel) = 0._dp
            moments_i(ip,ijeg_i,ikx,iky,iz,updatetlevel) = 0._dp
          ENDDO
          DO ij = ijsg_i,ijeg_i
            moments_i(ipsg_i+1,ij,ikx,iky,iz,updatetlevel) = 0._dp
            moments_i(ipsg_i  ,ij,ikx,iky,iz,updatetlevel) = 0._dp
            moments_i(ipeg_i-1,ij,ikx,iky,iz,updatetlevel) = 0._dp
            moments_i(ipeg_i  ,ij,ikx,iky,iz,updatetlevel) = 0._dp
          ENDDO
          kernel_i(ijsg_i,ikx,iky)      = 0._dp
          kernel_i(ijeg_i,ikx,iky)      = 0._dp
        ENDDO
      ENDDO
    ENDDO
  ELSE
    if(my_id .EQ. 0) write(*,*) '! Closure scheme not found !'

  ENDIF

  CALL cpu_time(t1_clos)
  tc_clos = tc_clos + (t1_clos - t0_clos)
END SUBROUTINE apply_closure_model

END module closure
