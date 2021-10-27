MODULE numerics
    USE basic
    USE prec_const
    USE grid
    USE utility
    USE coeff
    implicit none

    PUBLIC :: compute_derivatives, build_dnjs_table, evaluate_kernels, compute_lin_coeff
    PUBLIC :: wipe_turbulence, wipe_zonalflow

CONTAINS

! Compute the 2D particle temperature for electron and ions (sum over Laguerre)
SUBROUTINE compute_derivatives

END SUBROUTINE compute_derivatives

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
  USE array, Only : kernel_e, kernel_i, kparray
  USE grid
  USE model, ONLY : tau_e, tau_i, sigma_e, sigma_i, q_e, q_i, lambdaD, CLOS, sigmae2_taue_o2, sigmai2_taui_o2
  IMPLICIT NONE
  INTEGER    :: j_int
  REAL(dp)   :: j_dp, y_, kp2_, kx_, ky_

DO ikx = ikxs,ikxe
DO iky = ikys,ikye
DO iz = izs,ize
  !!!!! Electron kernels !!!!!
  DO ij = ijsg_e, ijeg_e
    j_int = jarray_e(ij)
    j_dp  = REAL(j_int,dp)
    y_    =  sigmae2_taue_o2 * kparray(ikx,iky,iz)**2
    kernel_e(ij,ikx,iky,iz) = y_**j_int*EXP(-y_)/GAMMA(j_dp+1._dp)!factj
  ENDDO
  !!!!! Ion kernels !!!!!
  DO ij = ijsg_i, ijeg_i
    j_int = jarray_i(ij)
    j_dp  = REAL(j_int,dp)
    y_    =  sigmai2_taui_o2 * kparray(ikx,iky,iz)**2
    kernel_i(ij,ikx,iky,iz) = y_**j_int*EXP(-y_)/GAMMA(j_dp+1._dp)!factj
  ENDDO
ENDDO
ENDDO
ENDDO

END SUBROUTINE evaluate_kernels
!******************************************************************************!

SUBROUTINE compute_lin_coeff
  USE array
  USE model, ONLY: taue_qe, taui_qi, sqrtTaue_qe, sqrtTaui_qi, &
                   K_T, K_n, CurvB, GradB
  USE prec_const
  USE grid,  ONLY: parray_e, parray_i, jarray_e, jarray_i, &
                   ip,ij, ips_e,ipe_e, ips_i,ipe_i, ijs_e,ije_e, ijs_i,ije_i
  IMPLICIT NONE
  INTEGER     :: p_int, j_int ! polynom. degrees
  REAL(dp)    :: p_dp, j_dp
  REAL(dp)    :: kx, ky, z
  !! Electrons linear coefficients for moment RHS !!!!!!!!!!
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
!!!!!!! Remove all ky==0 modes to conserve only non zonal modes in a restart
!******************************************************************************!
SUBROUTINE wipe_zonalflow
  USE fields
  USE grid
  IMPLICIT NONE
  DO ikx=ikxs,ikxe
  DO iky=ikys,ikye
  DO iz=izs,ize
    DO ip=ips_e,ipe_e
    DO ij=ijs_e,ije_e
      IF( iky .EQ. iky_0) THEN
        moments_e( ip,ij,ikx,iky,iz, :) = 0e-3_dp*moments_e( ip,ij,ikx,iky,iz, :)
      ELSE
        moments_e( ip,ij,ikx,iky,iz, :) = 1e+0_dp*moments_e( ip,ij,ikx,iky,iz, :)
      ENDIF
    ENDDO
    ENDDO
    DO ip=ips_i,ipe_i
    DO ij=ijs_i,ije_i
      IF( iky .EQ. iky_0) THEN
        moments_i( ip,ij,ikx,iky,iz, :) = 0e-3_dp*moments_i( ip,ij,ikx,iky,iz, :)
      ELSE
        moments_i( ip,ij,ikx,iky,iz, :) = 1e+0_dp*moments_i( ip,ij,ikx,iky,iz, :)
      ENDIF
    ENDDO
    ENDDO
    IF( iky .EQ. iky_0) THEN
      phi(ikx,iky,iz) = 0e-3_dp*phi(ikx,iky,iz)
    ELSE
      phi(ikx,iky,iz) = 1e+0_dp*phi(ikx,iky,iz)
    ENDIF
  ENDDO
  ENDDO
  ENDDO
END SUBROUTINE

END MODULE numerics
