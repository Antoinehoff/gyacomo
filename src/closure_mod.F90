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

SUBROUTINE apply_closure_model
  IMPLICIT NONE
  complex(dp) :: i_kz
  real(dp)    :: taue_qe_etaB_nu, taui_qi_etaB_nu
  real(dp)    :: sqpp2pp1_e, sqpp2pp1_i, sqpp1p_e, sqpp1p_i
  real(dp)    :: p_dp, j_dp
  ! Spare some computations
  taue_qe_etaB_nu = tau_e*eta_B/q_e/nu
  taui_qi_etaB_nu = tau_i*eta_B/q_i/nu
  sqpp2pp1_e = SQRT((pmaxe_dp+2)*(pmaxe_dp+1))
  sqpp2pp1_i = SQRT((pmaxi_dp+2)*(pmaxi_dp+1))
  sqpp1p_e   = SQRT((pmaxe_dp+1)*(pmaxe_dp))
  sqpp1p_i   = SQRT((pmaxi_dp+1)*(pmaxi_dp))

  CALL cpu_time(t0_clos)
    ! Positive Oob indices are approximated with a model
    IF (CLOS .EQ. 0) THEN
      ! zero truncation, An+1=0 for n+1>nmax
      DO ikr = ikrs,ikre
        DO ikz = ikzs,ikze

          DO ip = ipsg_e,ipeg_e
            moments_e(ip,ijsg_e,ikr,ikz,updatetlevel) = 0._dp
            moments_e(ip,ijeg_e,ikr,ikz,updatetlevel) = 0._dp
          ENDDO
          DO ij = ijsg_e,ijeg_e
            moments_e(ipsg_e+1,ij,ikr,ikz,updatetlevel) = 0._dp
            moments_e(ipsg_e  ,ij,ikr,ikz,updatetlevel) = 0._dp
            moments_e(ipeg_e-1,ij,ikr,ikz,updatetlevel) = 0._dp
            moments_e(ipeg_e  ,ij,ikr,ikz,updatetlevel) = 0._dp
          ENDDO
          kernel_e(ijsg_e,ikr,ikz)      = 0._dp
          kernel_e(ijeg_e,ikr,ikz)      = 0._dp

          DO ip = ipsg_i,ipeg_i
            moments_i(ip,ijsg_i,ikr,ikz,updatetlevel) = 0._dp
            moments_i(ip,ijeg_i,ikr,ikz,updatetlevel) = 0._dp
          ENDDO
          DO ij = ijsg_i,ijeg_i
            moments_i(ipsg_i+1,ij,ikr,ikz,updatetlevel) = 0._dp
            moments_i(ipsg_i  ,ij,ikr,ikz,updatetlevel) = 0._dp
            moments_i(ipeg_i-1,ij,ikr,ikz,updatetlevel) = 0._dp
            moments_i(ipeg_i  ,ij,ikr,ikz,updatetlevel) = 0._dp
          ENDDO
          kernel_i(ijsg_i,ikr,ikz)      = 0._dp
          kernel_i(ijeg_i,ikr,ikz)      = 0._dp

        ENDDO
      ENDDO
    ELSE
      if(my_id .EQ. 0) write(*,*) '! Closure scheme not found !'

    ENDIF

    CALL cpu_time(t1_clos)
    tc_clos = tc_clos + (t1_clos - t0_clos)
  END SUBROUTINE apply_closure_model
END module closure
