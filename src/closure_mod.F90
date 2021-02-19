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

  ! Negative out of bounds indices are put to zero (analytically correct)
    DO ikr = ikrs,ikre
      DO ikz = ikzs,ikze

        DO ip = ipsg_e,ipeg_e
          moments_e(ip,ijsg_e,ikr,ikz,:) = 0._dp
        ENDDO
        DO ij = ijsg_e,ijeg_e
          moments_e(ipsg_e+1,ij,ikr,ikz,:) = 0._dp
          moments_e(ipsg_e  ,ij,ikr,ikz,:) = 0._dp
        ENDDO
        kernel_e(ijsg_e,ikr,ikz)      = 0._dp

        DO ip = ipsg_i,ipeg_i
          moments_i(ip,ijsg_i,ikr,ikz,:) = 0._dp
        ENDDO
        DO ij = ijsg_i,ijeg_i
          moments_i(ipsg_i+1,ij,ikr,ikz,:) = 0._dp
          moments_i(ipsg_i  ,ij,ikr,ikz,:) = 0._dp
        ENDDO
        kernel_i(ijsg_i,ikr,ikz)      = 0._dp

      ENDDO
    ENDDO

    ! Positive Oob indices are approximated with a model
    IF (CLOS .EQ. 0) THEN
      ! zero truncation, An+1=0 for n+1>nmax
      DO ikr = ikrs,ikre
        DO ikz = ikzs,ikze

          DO ip = ipsg_e,ipeg_e
            moments_e(ip,ijeg_e,ikr,ikz,:) = 0._dp
          ENDDO
          DO ij = ijsg_e,ijeg_e
            moments_e(ipeg_e-1,ij,ikr,ikz,:) = 0._dp
            moments_e(ipeg_e  ,ij,ikr,ikz,:) = 0._dp
          ENDDO
          kernel_e(ijeg_e,ikr,ikz)      = 0._dp

          DO ip = ipsg_i,ipeg_i
            moments_i(ip,ijeg_i,ikr,ikz,:) = 0._dp
          ENDDO
          DO ij = ijsg_i,ijeg_i
            moments_i(ipeg_i-1,ij,ikr,ikz,:) = 0._dp
            moments_i(ipeg_i  ,ij,ikr,ikz,:) = 0._dp
          ENDDO
          kernel_i(ijeg_i,ikr,ikz)      = 0._dp
          
        ENDDO
      ENDDO

    ELSEIF ((CLOS .EQ. 1) .AND. (nu .NE. 0)) THEN
      !Semi collisional closure, i.e. at high degree, -nu N_M+1 ~ X_lin_M * N_M
      DO ikz = ikzs,ikze
        i_kz = imagu*kzarray(ikz)
        !! ELECTRONS
        ! Hermite closures
        DO ij = ijsg_e,ijeg_e
          j_dp = real(jarray_e(ij),dp)
          DO ikr = ikrs,ikre
            ! For p = Pmax + 2
            moments_e(ipeg_e,ij,ikr,ikz,:) = &
            -taue_qe_etaB_nu * i_kz * sqpp2pp1_e/(2*(pmaxe_dp+2)+j_dp) &
            * moments_e(ipe_e,ij,ikr,ikz,:)
            ! For p = Pmax + 1
            moments_e(ipeg_e-1,ij,ikr,ikz,:) = &
            -taue_qe_etaB_nu * i_kz * sqpp1p_e/(2*(pmaxe_dp+1)+j_dp) &
            * moments_e(ipe_e-1,ij,ikr,ikz,:)
            ! Kernel closure
          ENDDO
        ENDDO
        ! Laguerre closure
        DO ip = ipsg_e,ipeg_e
          p_dp = real(parray_e(ip),dp)
            DO ikr = ikrs,ikre
              ! For j = Jmax + 1
              moments_e(ip,ijeg_e,ikr,ikz,:) = &
              +taue_qe_etaB_nu * i_kz * (jmaxe_dp+1)/(2*p_dp+jmaxe_dp+1) &
              * moments_e(ip,ije_e,ikr,ikz,:)
            ENDDO
        ENDDO

        !! IONS
        ! Hermite closures
        DO ij = ijsg_i,ijeg_i
          j_dp = real(jarray_i(ij),dp)
          DO ikr = ikrs,ikre
            ! For p = Pmax + 2
            moments_i(ipeg_i,ij,ikr,ikz,:) = &
            -taui_qi_etaB_nu * i_kz * sqpp2pp1_i/(2*(pmaxi_dp+2)+j_dp) &
            * moments_i(ipe_i,ij,ikr,ikz,:)
            ! For p = Pmax + 1
            moments_i(ipeg_i-1,ij,ikr,ikz,:) = &
            -taui_qi_etaB_nu * i_kz * sqpp1p_i/(2*(pmaxi_dp+1)+j_dp) &
            * moments_i(ipe_i-1,ij,ikr,ikz,:)
          ENDDO
        ENDDO
        ! Laguerre closure
        DO ip = ipsg_i,ipeg_i
          p_dp = real(parray_i(ip),dp)
          DO ikr = ikrs,ikre
            ! For j = Jmax + 1
            moments_i(ip,ijeg_i,ikr,ikz,:) = &
            +taui_qi_etaB_nu * i_kz * (jmaxi_dp+1)/(2*p_dp+jmaxi_dp+1) &
            * moments_i(ip,ije_i,ikr,ikz,:)
          ENDDO
        ENDDO
      ENDDO

    ELSEIF (CLOS .EQ. 2) THEN
        ! Copy closure : P+2 <- P, P+1 <- P-1, J+1 <- J
        DO ikr = ikrs,ikre
          DO ikz = ikzs,ikze
  
            DO ip = ipsg_e,ipeg_e
              ! J ghost is J+1, so we put moment J to J+1
              moments_e(ip,ijeg_e,ikr,ikz,:) = moments_e(ip,ije_e,ikr,ikz,:)
            ENDDO
            DO ij = ijsg_e,ijeg_e
              ! P ghosts are P+1 and P+2, P+1 <- P-1 and P+2 <- P
              moments_e(ipeg_e-1,ij,ikr,ikz,:) = moments_e(ipe_e-1,ij,ikr,ikz,:)
              moments_e(ipeg_e  ,ij,ikr,ikz,:) = moments_e(ipe_e  ,ij,ikr,ikz,:)
            ENDDO
            ! Same for ions
            DO ip = ipsg_i,ipeg_i
              moments_i(ip,ijeg_i,ikr,ikz,:) = moments_i(ip,ije_i,ikr,ikz,:)
            ENDDO
            DO ij = ijsg_i,ijeg_i
              moments_i(ipeg_i-1,ij,ikr,ikz,:) = moments_i(ipe_i-1,ij,ikr,ikz,:)
              moments_i(ipeg_i  ,ij,ikr,ikz,:) = moments_i(ipe_i  ,ij,ikr,ikz,:)
            ENDDO
            
          ENDDO
        ENDDO

    ELSE
      if(my_id .EQ. 0) write(*,*) '! Closure scheme not found !'

    ENDIF

  END SUBROUTINE apply_closure_model
END module closure
