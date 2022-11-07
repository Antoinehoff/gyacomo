!!  The numerical_experiments module contains routines to "play" with the fourier
! modes in order to understand mechanisms. These routines are not integrated in
! the main code anymore because they are not used. This file serves as an archive.
MODULE numerical_experiments
USE basic
USE prec_const
USE grid
USE utility

implicit none

PUBLIC :: play_with_modes, save_EM_ZF_modes

CONTAINS
!******************************************************************************!
!!!!!!! Routine that can artificially increase or wipe modes
!******************************************************************************!
SUBROUTINE save_EM_ZF_modes
  USE fields
  USE array, ONLY : moments_e_ZF, moments_i_ZF, phi_ZF, moments_e_EM,moments_i_EM,phi_EM
  USE grid
  USE time_integration, ONLY: updatetlevel
  USE model, ONLY: KIN_E
  IMPLICIT NONE
  ! Store Zonal and entropy modes
  IF(contains_ky0) THEN
  IF(KIN_E) &
    moments_e_ZF(ips_e:ipe_e,ijs_e:ije_e,ikxs:ikxe,izs:ize) = moments_e(ips_e:ipe_e,ijs_e:ije_e,iky_0,ikxs:ikxe,izs:ize,updatetlevel)
    moments_i_ZF(ips_i:ipe_i,ijs_i:ije_i,ikxs:ikxe,izs:ize) = moments_i(ips_i:ipe_i,ijs_i:ije_i,iky_0,ikxs:ikxe,izs:ize,updatetlevel)
    phi_ZF(ikxs:ikxe,izs:ize) = phi(iky_0,ikxs:ikxe,izs:ize)
  ELSE
    IF(KIN_E) &
    moments_e_ZF(ips_e:ipe_e,ijs_e:ije_e,ikxs:ikxe,izs:ize) = 0._dp
    moments_i_ZF(ips_i:ipe_i,ijs_i:ije_i,ikxs:ikxe,izs:ize) = 0._dp
    phi_ZF(ikxs:ikxe,izs:ize) = 0._dp
  ENDIF
  IF(contains_kx0) THEN
    IF(KIN_E) &
    moments_e_EM(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,izs:ize) = moments_e(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,ikx_0,izs:ize,updatetlevel)
    moments_i_EM(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,izs:ize) = moments_i(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,ikx_0,izs:ize,updatetlevel)
    phi_EM(ikys:ikye,izs:ize) = phi(ikys:ikye,ikx_0,izs:ize)
  ELSE
    IF(KIN_E) &
    moments_e_EM(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,izs:ize) = 0._dp
    moments_i_EM(ips_e:ipe_e,ijs_i:ije_i,ikys:ikye,izs:ize) = 0._dp
    phi_EM(ikys:ikye,izs:ize) = 0._dp
  ENDIF
END SUBROUTINE

SUBROUTINE play_with_modes
  USE fields
  USE array, ONLY : moments_e_ZF, moments_i_ZF, phi_ZF, moments_e_EM,moments_i_EM,phi_EM
  USE grid
  USE time_integration, ONLY: updatetlevel
  USE initial_par, ONLY: ACT_ON_MODES
  USE model, ONLY: KIN_E
  IMPLICIT NONE
  REAL(dp) :: AMP = 1.5_dp

  SELECT CASE(ACT_ON_MODES)
  CASE('wipe_zonal') ! Errase the zonal flow
    IF(KIN_E) &
    moments_e(ips_e:ipe_e,ijs_e:ije_e,iky_0,ikxs:ikxe,izs:ize,updatetlevel) = 0._dp
    moments_i(ips_i:ipe_i,ijs_i:ije_i,iky_0,ikxs:ikxe,izs:ize,updatetlevel) = 0._dp
    phi(iky_0,ikxs:ikxe,izs:ize) = 0._dp
  CASE('wipe_entropymode')
    IF(KIN_E) &
    moments_e(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,ikx_0,izs:ize,updatetlevel) = 0._dp
    moments_i(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,ikx_0,izs:ize,updatetlevel) = 0._dp
    phi(ikys:ikye,ikx_0,izs:ize) = 0._dp
  CASE('wipe_turbulence')
    DO ikx = ikxs,ikxe
      DO iky = ikys, ikye
        IF ( (ikx .NE. ikx_0) .AND. (iky .NE. iky_0) ) THEN
          IF(KIN_E) &
          moments_e(ips_e:ipe_e,ijs_e:ije_e,iky,ikx,izs:ize,updatetlevel) = 0._dp
          moments_i(ips_i:ipe_i,ijs_i:ije_i,iky,ikx,izs:ize,updatetlevel) = 0._dp
          phi(iky,ikx,izs:ize) = 0._dp
        ENDIF
      ENDDO
    ENDDO
  CASE('wipe_nonzonal')
    DO ikx = ikxs,ikxe
      DO iky = ikys, ikye
        IF ( (ikx .NE. ikx_0) ) THEN
          IF(KIN_E) &
          moments_e(ips_e:ipe_e,ijs_e:ije_e,iky,ikx,izs:ize,updatetlevel) = 0._dp
          moments_i(ips_i:ipe_i,ijs_i:ije_i,iky,ikx,izs:ize,updatetlevel) = 0._dp
          phi(iky,ikx,izs:ize) = 0._dp
        ENDIF
      ENDDO
    ENDDO
  CASE('freeze_zonal')
    IF(KIN_E) &
    moments_e(ips_e:ipe_e,ijs_e:ije_e,iky_0,ikxs:ikxe,izs:ize,updatetlevel) = moments_e_ZF(ips_e:ipe_e,ijs_e:ije_e,ikxs:ikxe,izs:ize)
    moments_i(ips_i:ipe_i,ijs_i:ije_i,iky_0,ikxs:ikxe,izs:ize,updatetlevel) = moments_i_ZF(ips_i:ipe_i,ijs_i:ije_i,ikxs:ikxe,izs:ize)
    phi(iky_0,ikxs:ikxe,izs:ize) = phi_ZF(ikxs:ikxe,izs:ize)
  CASE('freeze_entropymode')
    IF(contains_kx0) THEN
      IF(KIN_E) &
      moments_e(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,ikx_0,izs:ize,updatetlevel) = moments_e_EM(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,izs:ize)
      moments_i(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,ikx_0,izs:ize,updatetlevel) = moments_i_EM(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,izs:ize)
      phi(ikys:ikye,ikx_0,izs:ize) = phi_EM(ikys:ikye,izs:ize)
    ENDIF
  CASE('amplify_zonal')
    IF(KIN_E) &
    moments_e(ips_e:ipe_e,ijs_e:ije_e,iky_0,ikxs:ikxe,izs:ize,updatetlevel) = AMP*moments_e_ZF(ips_e:ipe_e,ijs_e:ije_e,ikxs:ikxe,izs:ize)
    moments_i(ips_i:ipe_i,ijs_i:ije_i,iky_0,ikxs:ikxe,izs:ize,updatetlevel) = AMP*moments_i_ZF(ips_i:ipe_i,ijs_i:ije_i,ikxs:ikxe,izs:ize)
    phi(iky_0,ikxs:ikxe,izs:ize) = AMP*phi_ZF(ikxs:ikxe,izs:ize)
  END SELECT
END SUBROUTINE

END MODULE numerical experiments
