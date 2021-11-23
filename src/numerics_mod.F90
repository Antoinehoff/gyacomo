MODULE numerics
    USE basic
    USE prec_const
    USE grid
    USE utility
    USE coeff
    implicit none

    PUBLIC :: build_dnjs_table, evaluate_kernels, evaluate_poisson_op, compute_lin_coeff
    PUBLIC :: wipe_turbulence, play_with_modes, save_EM_ZF_modes

CONTAINS

!******************************************************************************!
!!!!!!! Build the Laguerre-Laguerre coupling coefficient table for nonlin
!******************************************************************************!
SUBROUTINE build_dnjs_table
  USE array, Only : dnjs
  USE coeff
  IMPLICIT NONE

  INTEGER :: in, ij, is, J
  INTEGER :: n_, j_, s_

  J = max(jmaxe,jmaxi)

  DO in = 1,J+1 ! Nested dependent loops to make benefit from dnjs symmetry
    n_ = in - 1
    DO ij = in,J+1
      j_ = ij - 1
      DO is = ij,J+1
        s_ = is - 1

        dnjs(in,ij,is) = TO_DP(ALL2L(n_,j_,s_,0))
        ! By symmetry
        dnjs(in,is,ij) = dnjs(in,ij,is)
        dnjs(ij,in,is) = dnjs(in,ij,is)
        dnjs(ij,is,in) = dnjs(in,ij,is)
        dnjs(is,ij,in) = dnjs(in,ij,is)
        dnjs(is,in,ij) = dnjs(in,ij,is)
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE build_dnjs_table
!******************************************************************************!

!******************************************************************************!
!!!!!!! Evaluate the kernels once for all
!******************************************************************************!
SUBROUTINE evaluate_kernels
  USE basic
  USE array, Only : kernel_e, kernel_i
  USE grid
  USE model, ONLY : tau_e, tau_i, sigma_e, sigma_i, q_e, q_i, &
                    lambdaD, CLOS, sigmae2_taue_o2, sigmai2_taui_o2, KIN_E
  IMPLICIT NONE
  INTEGER    :: j_int
  REAL(dp)   :: j_dp, y_, kp2_, kx_, ky_

DO eo  = 0,1
DO ikx = ikxs,ikxe
DO iky = ikys,ikye
DO iz = izs,ize
  !!!!! Electron kernels !!!!!
  IF(KIN_E) THEN
  DO ij = ijsg_e, ijeg_e
    j_int = jarray_e(ij)
    j_dp  = REAL(j_int,dp)
    y_    =  sigmae2_taue_o2 * kparray(ikx,iky,iz,eo)**2
    kernel_e(ij,ikx,iky,iz,eo) = y_**j_int*EXP(-y_)/GAMMA(j_dp+1._dp)!factj
  ENDDO
  ENDIF
  !!!!! Ion kernels !!!!!
  DO ij = ijsg_i, ijeg_i
    j_int = jarray_i(ij)
    j_dp  = REAL(j_int,dp)
    y_    =  sigmai2_taui_o2 * kparray(ikx,iky,iz,eo)**2
    kernel_i(ij,ikx,iky,iz,eo) = y_**j_int*EXP(-y_)/GAMMA(j_dp+1._dp)!factj
  ENDDO
ENDDO
ENDDO
ENDDO
ENDDO

END SUBROUTINE evaluate_kernels
!******************************************************************************!

!******************************************************************************!
!!!!!!! Evaluate polarisation operator for Poisson equation
!******************************************************************************!
SUBROUTINE evaluate_poisson_op
  USE basic
  USE array, Only : kernel_e, kernel_i, inv_poisson_op
  USE grid
  USE model, ONLY : tau_e, tau_i, q_e, q_i, KIN_E
  IMPLICIT NONE
  REAL(dp)    :: pol_i, pol_e     ! (Z_a^2/tau_a (1-sum_n kernel_na^2))
  INTEGER     :: ini,ine

  kxloop: DO ikx = ikxs,ikxe
  kyloop: DO iky = ikys,ikye
  zloop:  DO iz  =  izs,ize
  ! This term is evalued on the even z grid since poisson uses only p=0 and phi
  IF( (kxarray(ikx).EQ.0._dp) .AND. (kyarray(iky).EQ.0._dp) ) THEN
      inv_poisson_op(ikx, iky, iz) =  0._dp
    ELSE
    !!!!!!!!!!!!!!!!! Ion contribution
    ! loop over n only if the max polynomial degree
    pol_i = 0._dp
    DO ini=1,jmaxi+1
      pol_i = pol_i  + qi2_taui*kernel_i(ini,ikx,iky,iz,0)**2 ! ... sum recursively ...
    END DO
    !!!!!!!!!!!!! Electron contribution\
    pol_e = 0._dp
    ! Kinetic model
    IF (KIN_E) THEN
      ! loop over n only if the max polynomial degree
      DO ine=1,jmaxe+1 ! ine = n+1
        pol_e = pol_e  + qe2_taue*kernel_e(ine,ikx,iky,iz,0)**2 ! ... sum recursively ...
      END DO
    ! Adiabatic model
    ELSE
      pol_e = 1._dp - qe2_taue
    ENDIF
    inv_poisson_op(ikx, iky, iz) =  1._dp/(qe2_taue + qi2_taui - pol_i - pol_e)
  ENDIF
  END DO zloop
  END DO kyloop
  END DO kxloop

END SUBROUTINE evaluate_poisson_op
!******************************************************************************!

SUBROUTINE compute_lin_coeff
  USE array
  USE model, ONLY: taue_qe, taui_qi, sqrtTaue_qe, sqrtTaui_qi, &
                   K_T, K_n, CurvB, GradB, KIN_E
  USE prec_const
  USE grid,  ONLY: parray_e, parray_i, jarray_e, jarray_i, &
                   ip,ij, ips_e,ipe_e, ips_i,ipe_i, ijs_e,ije_e, ijs_i,ije_i
  IMPLICIT NONE
  INTEGER     :: p_int, j_int ! polynom. degrees
  REAL(dp)    :: p_dp, j_dp
  REAL(dp)    :: kx, ky, z
  !! Electrons linear coefficients for moment RHS !!!!!!!!!!
  IF(KIN_E)THEN
  DO ip = ips_e, ipe_e
    p_int= parray_e(ip)   ! Hermite degree
    p_dp = REAL(p_int,dp) ! REAL of Hermite degree
    DO ij = ijs_e, ije_e
      j_int= jarray_e(ij)   ! Laguerre degree
      j_dp = REAL(j_int,dp) ! REAL of Laguerre degree
      ! All Napj terms
      xnepj(ip,ij) = taue_qe*(CurvB*(2._dp*p_dp + 1._dp) &
                             +GradB*(2._dp*j_dp + 1._dp))
      ! Mirror force terms
      ynepp1j  (ip,ij) = -SQRT(tau_e)/sigma_e *          (j_dp+1)*SQRT(p_dp+1._dp)
      ynepm1j  (ip,ij) = -SQRT(tau_e)/sigma_e *          (j_dp+1)*SQRT(p_dp)
      ynepp1jm1(ip,ij) = +SQRT(tau_e)/sigma_e *              j_dp*SQRT(p_dp+1._dp)
      ynepm1jm1(ip,ij) = +SQRT(tau_e)/sigma_e *              j_dp*SQRT(p_dp)
      zNepm1j  (ip,ij) = +SQRT(tau_e)/sigma_e * (2._dp*j_dp+1_dp)*SQRT(p_dp)
      zNepm1jp1(ip,ij) = -SQRT(tau_e)/sigma_e *       (j_dp+1_dp)*SQRT(p_dp)
      zNepm1jm1(ip,ij) = -SQRT(tau_e)/sigma_e *              j_dp*SQRT(p_dp)
    ENDDO
  ENDDO
  DO ip = ips_e, ipe_e
    p_int= parray_e(ip)   ! Hermite degree
    p_dp = REAL(p_int,dp) ! REAL of Hermite degree
    ! Landau damping coefficients (ddz napj term)
    xnepp1j(ip) = SQRT(tau_e)/sigma_e * SQRT(p_dp + 1_dp)
    xnepm1j(ip) = SQRT(tau_e)/sigma_e * SQRT(p_dp)
    ! Magnetic curvature term
    xnepp2j(ip) = taue_qe * CurvB * SQRT((p_dp + 1._dp) * (p_dp + 2._dp))
    xnepm2j(ip) = taue_qe * CurvB * SQRT(p_dp * (p_dp - 1._dp))
  ENDDO
  DO ij = ijs_e, ije_e
    j_int= jarray_e(ij)   ! Laguerre degree
    j_dp = REAL(j_int,dp) ! REAL of Laguerre degree
    ! Magnetic gradient term
    xnepjp1(ij) = -taue_qe * GradB * (j_dp + 1._dp)
    xnepjm1(ij) = -taue_qe * GradB * j_dp
  ENDDO
  ENDIF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Ions linear coefficients for moment RHS !!!!!!!!!!
  DO ip = ips_i, ipe_i
    p_int= parray_i(ip)   ! Hermite degree
    p_dp = REAL(p_int,dp) ! REAL of Hermite degree
    DO ij = ijs_i, ije_i
      j_int= jarray_i(ij)   ! Laguerre degree
      j_dp = REAL(j_int,dp) ! REAL of Laguerre degree
      ! All Napj terms
      xnipj(ip,ij) = taui_qi*(CurvB*(2._dp*p_dp + 1._dp) &
                             +GradB*(2._dp*j_dp + 1._dp))
      ! Mirror force terms
      ynipp1j  (ip,ij) = -SQRT(tau_i)/sigma_i*          (j_dp+1)*SQRT(p_dp+1._dp)
      ynipm1j  (ip,ij) = -SQRT(tau_i)/sigma_i*          (j_dp+1)*SQRT(p_dp)
      ynipp1jm1(ip,ij) = +SQRT(tau_i)/sigma_i*              j_dp*SQRT(p_dp+1._dp)
      ynipm1jm1(ip,ij) = +SQRT(tau_i)/sigma_i*              j_dp*SQRT(p_dp)
      ! Trapping terms
      zNipm1j  (ip,ij) = +SQRT(tau_i)/sigma_i* (2._dp*j_dp+1_dp)*SQRT(p_dp)
      zNipm1jp1(ip,ij) = -SQRT(tau_i)/sigma_i*       (j_dp+1_dp)*SQRT(p_dp)
      zNipm1jm1(ip,ij) = -SQRT(tau_i)/sigma_i*              j_dp*SQRT(p_dp)
    ENDDO
  ENDDO
  DO ip = ips_i, ipe_i
    p_int= parray_i(ip)   ! Hermite degree
    p_dp = REAL(p_int,dp) ! REAL of Hermite degree
    ! Landau damping coefficients (ddz napj term)
    xnipp1j(ip) = SQRT(tau_i)/sigma_i * SQRT(p_dp + 1._dp)
    xnipm1j(ip) = SQRT(tau_i)/sigma_i * SQRT(p_dp)
    ! Magnetic curvature term
    xnipp2j(ip) = taui_qi * CurvB * SQRT((p_dp + 1._dp) * (p_dp + 2._dp))
    xnipm2j(ip) = taui_qi * CurvB * SQRT(p_dp * (p_dp - 1._dp))
  ENDDO
  DO ij = ijs_i, ije_i
    j_int= jarray_i(ij)   ! Laguerre degree
    j_dp = REAL(j_int,dp) ! REAL of Laguerre degree
    ! Magnetic gradient term
    xnipjp1(ij) = -taui_qi * GradB * (j_dp + 1._dp)
    xnipjm1(ij) = -taui_qi * GradB * j_dp
  ENDDO
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! ES linear coefficients for moment RHS !!!!!!!!!!
  DO ip = ips_i, ipe_i
    p_int= parray_i(ip)   ! Hermite degree
    DO ij = ijs_i, ije_i
      j_int= jarray_i(ij)   ! REALof Laguerre degree
      j_dp = REAL(j_int,dp) ! REALof Laguerre degree
      !! Electrostatic potential pj terms
      IF (p_int .EQ. 0) THEN ! kronecker p0
        xphij(ip,ij)    =+K_n + 2.*j_dp*K_T
        xphijp1(ip,ij)  =-K_T*(j_dp+1._dp)
        xphijm1(ip,ij)  =-K_T* j_dp
      ELSE IF (p_int .EQ. 2) THEN ! kronecker p2
        xphij(ip,ij)    =+K_T/SQRT2
        xphijp1(ip,ij)  = 0._dp; xphijm1(ip,ij)  = 0._dp;
      ELSE
        xphij(ip,ij)    = 0._dp; xphijp1(ip,ij)  = 0._dp
        xphijm1(ip,ij)  = 0._dp;
      ENDIF
    ENDDO
  ENDDO

END SUBROUTINE compute_lin_coeff

!******************************************************************************!
!!!!!!! Remove all ky!=0 modes to conserve only zonal modes in a restart
!******************************************************************************!
SUBROUTINE wipe_turbulence
  USE fields
  USE grid
  IMPLICIT NONE
  DO ikx=ikxs,ikxe
  DO iky=ikys,ikye
  DO iz=izs,ize
    DO ip=ips_e,ipe_e
    DO ij=ijs_e,ije_e
      IF( iky .NE. iky_0) THEN
        moments_e( ip,ij,ikx,iky,iz, :) = 0e-3_dp*moments_e( ip,ij,ikx,iky,iz, :)
      ELSE
        moments_e( ip,ij,ikx,iky,iz, :) = 1e+0_dp*moments_e( ip,ij,ikx,iky,iz, :)
      ENDIF
    ENDDO
    ENDDO
    DO ip=ips_i,ipe_i
    DO ij=ijs_i,ije_i
      IF( iky .NE. iky_0) THEN
        moments_i( ip,ij,ikx,iky,iz, :) = 0e-3_dp*moments_i( ip,ij,ikx,iky,iz, :)
      ELSE
        moments_i( ip,ij,ikx,iky,iz, :) = 1e+0_dp*moments_i( ip,ij,ikx,iky,iz, :)
      ENDIF
    ENDDO
    ENDDO
    IF( iky .NE. iky_0) THEN
      phi(ikx,iky,iz) = 0e-3_dp*phi(ikx,iky,iz)
    ELSE
      phi(ikx,iky,iz) = 1e+0_dp*phi(ikx,iky,iz)
    ENDIF
  ENDDO
  ENDDO
  ENDDO
END SUBROUTINE
!******************************************************************************!
!!!!!!! Routine that can artificially increase or wipe modes
!******************************************************************************!
SUBROUTINE save_EM_ZF_modes
  USE fields
  USE array, ONLY : moments_e_ZF, moments_i_ZF, phi_ZF, moments_e_EM,moments_i_EM,phi_EM
  USE grid
  USE time_integration, ONLY: updatetlevel
  USE initial_par, ONLY: ACT_ON_MODES
  IMPLICIT NONE
  ! Store Zonal and entropy modes
  moments_e_ZF(ips_e:ipe_e,ijs_e:ije_e,ikxs:ikxe,izs:ize) = moments_e(ips_e:ipe_e,ijs_e:ije_e,ikxs:ikxe,iky_0,izs:ize,updatetlevel)
  moments_i_ZF(ips_i:ipe_i,ijs_i:ije_i,ikxs:ikxe,izs:ize) = moments_i(ips_i:ipe_i,ijs_i:ije_i,ikxs:ikxe,iky_0,izs:ize,updatetlevel)
  phi_ZF(ikxs:ikxe,izs:ize) = phi(ikxs:ikxe,iky_0,izs:ize)
  IF(contains_kx0) THEN
    moments_e_EM(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,izs:ize) = moments_e(ips_e:ipe_e,ijs_e:ije_e,ikx_0,ikys:ikye,izs:ize,updatetlevel)
    moments_i_EM(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,izs:ize) = moments_i(ips_i:ipe_i,ijs_i:ije_i,ikx_0,ikys:ikye,izs:ize,updatetlevel)
    phi_EM(ikys:ikye,izs:ize) = phi(ikx_0,ikys:ikye,izs:ize)
  ELSE
    moments_e_EM(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,izs:ize) = 0._dp
    moments_i_EM(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,izs:ize) = 0._dp
    phi_EM(ikys:ikye,izs:ize) = 0._dp
  ENDIF
END SUBROUTINE

SUBROUTINE play_with_modes
  USE fields
  USE array, ONLY : moments_e_ZF, moments_i_ZF, phi_ZF, moments_e_EM,moments_i_EM,phi_EM
  USE grid
  USE time_integration, ONLY: updatetlevel
  USE initial_par, ONLY: ACT_ON_MODES
  IMPLICIT NONE
  REAL(dp) :: AMP = 1.5_dp

  SELECT CASE(ACT_ON_MODES)
  CASE('wipe_zonal') ! Errase the zonal flow
    moments_e(ips_e:ipe_e,ijs_e:ije_e,ikxs:ikxe,iky_0,izs:ize,updatetlevel) = 0._dp
    moments_i(ips_i:ipe_i,ijs_i:ije_i,ikxs:ikxe,iky_0,izs:ize,updatetlevel) = 0._dp
    phi(ikxs:ikxe,iky_0,izs:ize) = 0._dp
  CASE('wipe_entropymode')
    moments_e(ips_e:ipe_e,ijs_e:ije_e,ikx_0,ikys:ikye,izs:ize,updatetlevel) = 0._dp
    moments_i(ips_i:ipe_i,ijs_i:ije_i,ikx_0,ikys:ikye,izs:ize,updatetlevel) = 0._dp
    phi(ikx_0,ikys:ikye,izs:ize) = 0._dp
  CASE('wipe_turbulence')
    DO ikx = ikxs,ikxe
      DO iky = ikys, ikye
        IF ( (ikx .NE. ikx_0) .AND. (iky .NE. iky_0) ) THEN
          moments_e(ips_e:ipe_e,ijs_e:ije_e,ikx,iky,izs:ize,updatetlevel) = 0._dp
          moments_i(ips_i:ipe_i,ijs_i:ije_i,ikx,iky,izs:ize,updatetlevel) = 0._dp
          phi(ikx,iky,izs:ize) = 0._dp
        ENDIF
      ENDDO
    ENDDO
  CASE('wipe_nonzonal')
    DO ikx = ikxs,ikxe
      DO iky = ikys, ikye
        IF ( (ikx .NE. ikx_0) ) THEN
          moments_e(ips_e:ipe_e,ijs_e:ije_e,ikx,iky,izs:ize,updatetlevel) = 0._dp
          moments_i(ips_i:ipe_i,ijs_i:ije_i,ikx,iky,izs:ize,updatetlevel) = 0._dp
          phi(ikx,iky,izs:ize) = 0._dp
        ENDIF
      ENDDO
    ENDDO
  CASE('freeze_zonal')
    moments_e(ips_e:ipe_e,ijs_e:ije_e,ikxs:ikxe,iky_0,izs:ize,updatetlevel) = moments_e_ZF(ips_e:ipe_e,ijs_e:ije_e,ikxs:ikxe,izs:ize)
    moments_i(ips_i:ipe_i,ijs_i:ije_i,ikxs:ikxe,iky_0,izs:ize,updatetlevel) = moments_i_ZF(ips_i:ipe_i,ijs_i:ije_i,ikxs:ikxe,izs:ize)
    phi(ikxs:ikxe,iky_0,izs:ize) = phi_ZF(ikxs:ikxe,izs:ize)
  CASE('freeze_entropymode')
    IF(contains_kx0) THEN
      moments_e(ips_e:ipe_e,ijs_e:ije_e,ikx_0,ikys:ikye,izs:ize,updatetlevel) = moments_e_EM(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,izs:ize)
      moments_i(ips_i:ipe_i,ijs_i:ije_i,ikx_0,ikys:ikye,izs:ize,updatetlevel) = moments_i_EM(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,izs:ize)
      phi(ikx_0,ikys:ikye,izs:ize) = phi_EM(ikys:ikye,izs:ize)
    ENDIF
  CASE('amplify_zonal')
    moments_e(ips_e:ipe_e,ijs_e:ije_e,ikxs:ikxe,iky_0,izs:ize,updatetlevel) = AMP*moments_e_ZF(ips_e:ipe_e,ijs_e:ije_e,ikxs:ikxe,izs:ize)
    moments_i(ips_i:ipe_i,ijs_i:ije_i,ikxs:ikxe,iky_0,izs:ize,updatetlevel) = AMP*moments_i_ZF(ips_i:ipe_i,ijs_i:ije_i,ikxs:ikxe,izs:ize)
    phi(ikxs:ikxe,iky_0,izs:ize) = AMP*phi_ZF(ikxs:ikxe,izs:ize)
  END SELECT
END SUBROUTINE

END MODULE numerics
