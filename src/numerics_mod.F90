MODULE numerics
    USE basic
    USE prec_const
    USE grid
    USE utility
    USE coeff
    implicit none

    PUBLIC :: build_dnjs_table, evaluate_kernels, evaluate_poisson_op, evaluate_ampere_op
    PUBLIC :: compute_lin_coeff, play_with_modes, save_EM_ZF_modes

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
  USE array, Only : kernel_e, kernel_i, HF_phi_correction_operator
  USE grid
  USE model, ONLY : tau_e, tau_i, sigma_e, sigma_i, q_e, q_i, &
                    lambdaD, CLOS, sigmae2_taue_o2, sigmai2_taui_o2, KIN_E
  IMPLICIT NONE
  INTEGER    :: j_int
  REAL(dp)   :: j_dp, y_

DO eo  = 0,1
DO ikx = ikxs,ikxe
DO iky = ikys,ikye
DO iz  = izgs,izge
  !!!!! Electron kernels !!!!!
  IF(KIN_E) THEN
  DO ij = ijgs_e, ijge_e
    j_int = jarray_e(ij)
    j_dp  = REAL(j_int,dp)
    y_    =  sigmae2_taue_o2 * kparray(iky,ikx,iz,eo)**2
    kernel_e(ij,iky,ikx,iz,eo) = y_**j_int*EXP(-y_)/GAMMA(j_dp+1._dp)!factj
  ENDDO
  IF (ijs_e .EQ. 1) &
  kernel_e(ijgs_e,iky,ikx,iz,eo) = 0._dp
  ENDIF
  !!!!! Ion kernels !!!!!
  DO ij = ijgs_i, ijge_i
    j_int = jarray_i(ij)
    j_dp  = REAL(j_int,dp)
    y_    =  sigmai2_taui_o2 * kparray(iky,ikx,iz,eo)**2
    kernel_i(ij,iky,ikx,iz,eo) = y_**j_int*EXP(-y_)/GAMMA(j_dp+1._dp)!factj
  ENDDO
  IF (ijs_i .EQ. 1) &
  kernel_i(ijgs_i,iky,ikx,iz,eo) = 0._dp
ENDDO
ENDDO
ENDDO
ENDDO
!! Correction term for the evaluation of the heat flux
HF_phi_correction_operator(ikys:ikye,ikxs:ikxe,izs:ize) = &
       2._dp * Kernel_i(1,ikys:ikye,ikxs:ikxe,izs:ize,0) &
      -1._dp * Kernel_i(2,ikys:ikye,ikxs:ikxe,izs:ize,0)

DO ij = ijs_i, ije_i
  j_int = jarray_i(ij)
  j_dp  = REAL(j_int,dp)
  HF_phi_correction_operator(ikys:ikye,ikxs:ikxe,izs:ize) = HF_phi_correction_operator(ikys:ikye,ikxs:ikxe,izs:ize) &
  - Kernel_i(ij,ikys:ikye,ikxs:ikxe,izs:ize,0) * (&
      2._dp*(j_dp+1.5_dp) * Kernel_i(ij  ,ikys:ikye,ikxs:ikxe,izs:ize,0) &
      -     (j_dp+1.0_dp) * Kernel_i(ij+1,ikys:ikye,ikxs:ikxe,izs:ize,0) &
      -              j_dp * Kernel_i(ij-1,ikys:ikye,ikxs:ikxe,izs:ize,0))
ENDDO

END SUBROUTINE evaluate_kernels
!******************************************************************************!

!******************************************************************************!
!!!!!!! Evaluate inverse polarisation operator for Poisson equation
!******************************************************************************!
SUBROUTINE evaluate_poisson_op
  USE basic
  USE array, Only : kernel_e, kernel_i, inv_poisson_op, inv_pol_ion
  USE grid
  USE model, ONLY : tau_e, tau_i, q_e, q_i, KIN_E
  IMPLICIT NONE
  REAL(dp)    :: pol_i, pol_e     ! (Z_a^2/tau_a (1-sum_n kernel_na^2))
  INTEGER     :: ini,ine

  ! This term has no staggered grid dependence. It is evalued for the
  ! even z grid since poisson uses p=0 moments and phi only.
  kxloop: DO ikx = ikxs,ikxe
  kyloop: DO iky = ikys,ikye
  zloop:  DO iz  = izs,ize
  IF( (kxarray(ikx).EQ.0._dp) .AND. (kyarray(iky).EQ.0._dp) ) THEN
      inv_poisson_op(iky, ikx, iz) =  0._dp
  ELSE
    !!!!!!!!!!!!!!!!! Ion contribution
    ! loop over n only if the max polynomial degree
    pol_i = 0._dp
    DO ini=1,jmaxi+1
      pol_i = pol_i  + qi2_taui*kernel_i(ini,iky,ikx,iz,0)**2 ! ... sum recursively ...
    END DO
    !!!!!!!!!!!!! Electron contribution
    pol_e = 0._dp
    ! Kinetic model
    IF (KIN_E) THEN
      ! loop over n only if the max polynomial degree
      DO ine=1,jmaxe+1 ! ine = n+1
        pol_e = pol_e  + qe2_taue*kernel_e(ine,iky,ikx,iz,0)**2 ! ... sum recursively ...
      END DO
    ! Adiabatic model
    ELSE
      pol_e = qe2_taue - 1._dp
    ENDIF
    inv_poisson_op(iky, ikx, iz) =  1._dp/(qe2_taue + qi2_taui - pol_i - pol_e)
    inv_pol_ion   (iky, ikx, iz) =  1._dp/(qi2_taui - pol_i)
  ENDIF
  END DO zloop
  END DO kyloop
  END DO kxloop

END SUBROUTINE evaluate_poisson_op
!******************************************************************************!

!******************************************************************************!
!!!!!!! Evaluate inverse polarisation operator for Poisson equation
!******************************************************************************!
SUBROUTINE evaluate_ampere_op
  USE basic
  USE array, Only : kernel_e, kernel_i, inv_ampere_op
  USE grid
  USE model, ONLY : tau_e, tau_i, q_e, q_i, KIN_E, beta
  IMPLICIT NONE
  REAL(dp)    :: pol_i, pol_e, kperp2     ! (Z_a^2/tau_a (1-sum_n kernel_na^2))
  INTEGER     :: ini,ine

  ! We do not solve Ampere if beta = 0 to spare waste of ressources
  IF(SOLVE_AMPERE) THEN
    ! This term has no staggered grid dependence. It is evalued for the
    ! even z grid since poisson uses p=0 moments and phi only.
    kxloop: DO ikx = ikxs,ikxe
    kyloop: DO iky = ikys,ikye
    zloop:  DO iz  = izs,ize
    kperp2 = kparray(iky,ikx,iz,0)**2
    IF( (kxarray(ikx).EQ.0._dp) .AND. (kyarray(iky).EQ.0._dp) ) THEN
        inv_ampere_op(iky, ikx, iz) =  0._dp
    ELSE
      !!!!!!!!!!!!!!!!! Ion contribution
      ! loop over n only if the max polynomial degree
      pol_i = 0._dp
      DO ini=1,jmaxi+1
        pol_i = pol_i  + kernel_i(ini,iky,ikx,iz,0)**2 ! ... sum recursively ...
      END DO
      pol_i = q_i**2/(sigma_i**2) * pol_i
      !!!!!!!!!!!!! Electron contribution
      pol_e = 0._dp
      ! loop over n only if the max polynomial degree
      DO ine=1,jmaxe+1 ! ine = n+1
        pol_e = pol_e  + kernel_e(ine,iky,ikx,iz,0)**2 ! ... sum recursively ...
      END DO
      pol_e = q_e**2/(sigma_e**2) * pol_e
      inv_ampere_op(iky, ikx, iz) =  1._dp/(2._dp*kperp2 + beta*(pol_i + pol_e))
    ENDIF
    END DO zloop
    END DO kyloop
    END DO kxloop
  ENDIF

END SUBROUTINE evaluate_ampere_op
!******************************************************************************!

SUBROUTINE compute_lin_coeff
  USE array
  USE model, ONLY: taue_qe, taui_qi, sqrtTaue_qe, sqrtTaui_qi, &
                   k_T, eta_T, k_N, eta_N, CurvB, GradB, KIN_E,&
                   tau_e, tau_i, sigma_e, sigma_i
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
  IF (KIN_E) THEN
    DO ip = ips_e, ipe_e
      p_int= parray_e(ip)   ! Hermite degree
      DO ij = ijs_e, ije_e
        j_int= jarray_e(ij)   ! REALof Laguerre degree
        j_dp = REAL(j_int,dp) ! REALof Laguerre degree
        !! Electrostatic potential pj terms
        IF (p_int .EQ. 0) THEN ! kronecker p0
          xphij_e(ip,ij)    =+eta_N*k_N + 2.*j_dp*eta_T*k_T
          xphijp1_e(ip,ij)  =-eta_T*k_T*(j_dp+1._dp)
          xphijm1_e(ip,ij)  =-eta_T*k_T* j_dp
        ELSE IF (p_int .EQ. 2) THEN ! kronecker p2
          xphij_e(ip,ij)    =+eta_T*k_T/SQRT2
          xphijp1_e(ip,ij)  = 0._dp; xphijm1_e(ip,ij)  = 0._dp;
        ELSE
          xphij_e(ip,ij)    = 0._dp; xphijp1_e(ip,ij)  = 0._dp
          xphijm1_e(ip,ij)  = 0._dp;
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  DO ip = ips_i, ipe_i
    p_int= parray_i(ip)   ! Hermite degree
    DO ij = ijs_i, ije_i
      j_int= jarray_i(ij)   ! REALof Laguerre degree
      j_dp = REAL(j_int,dp) ! REALof Laguerre degree
      !! Electrostatic potential pj terms
      IF (p_int .EQ. 0) THEN ! kronecker p0
        xphij_i(ip,ij)    =+k_N + 2._dp*j_dp*k_T
        xphijp1_i(ip,ij)  =-k_T*(j_dp+1._dp)
        xphijm1_i(ip,ij)  =-k_T* j_dp
      ELSE IF (p_int .EQ. 2) THEN ! kronecker p2
        xphij_i(ip,ij)    =+k_T/SQRT2
        xphijp1_i(ip,ij)  = 0._dp; xphijm1_i(ip,ij)  = 0._dp;
      ELSE
        xphij_i(ip,ij)    = 0._dp; xphijp1_i(ip,ij)  = 0._dp
        xphijm1_i(ip,ij)  = 0._dp;
      ENDIF
    ENDDO
  ENDDO
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! EM linear coefficients for moment RHS !!!!!!!!!!
  IF (KIN_E) THEN
    DO ip = ips_e, ipe_e
      p_int= parray_e(ip)   ! Hermite degree
      DO ij = ijs_e, ije_e
        j_int= jarray_e(ij)   ! REALof Laguerre degree
        j_dp = REAL(j_int,dp) ! REALof Laguerre degree
        !! Electrostatic potential pj terms
        IF (p_int .EQ. 1) THEN ! kronecker p1
          xpsij_e(ip,ij)    =+(eta_N*k_N + (2._dp*j_dp+1._dp)*eta_T*k_T) * SQRT(tau_e)/sigma_e
          xpsijp1_e(ip,ij)  =- eta_T*k_T*(j_dp+1._dp) * SQRT(tau_e)/sigma_e
          xpsijm1_e(ip,ij)  =- eta_T*k_T* j_dp        * SQRT(tau_e)/sigma_e
        ELSE IF (p_int .EQ. 3) THEN ! kronecker p3
          xpsij_e(ip,ij)    =+ eta_T*k_T*SQRT3/SQRT2  * SQRT(tau_e)/sigma_e
          xpsijp1_e(ip,ij)  = 0._dp; xpsijm1_e(ip,ij)  = 0._dp;
        ELSE
          xpsij_e(ip,ij)    = 0._dp; xpsijp1_e(ip,ij)  = 0._dp
          xpsijm1_e(ip,ij)  = 0._dp;
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  DO ip = ips_i, ipe_i
    p_int= parray_i(ip)   ! Hermite degree
    DO ij = ijs_i, ije_i
      j_int= jarray_i(ij)   ! REALof Laguerre degree
      j_dp = REAL(j_int,dp) ! REALof Laguerre degree
      !! Electrostatic potential pj terms
      IF (p_int .EQ. 1) THEN ! kronecker p1
        xpsij_i(ip,ij)    =+(k_N + (2._dp*j_dp+1._dp)*k_T) * SQRT(tau_i)/sigma_i
        xpsijp1_i(ip,ij)  =- k_T*(j_dp+1._dp)              * SQRT(tau_i)/sigma_i
        xpsijm1_i(ip,ij)  =- k_T* j_dp                     * SQRT(tau_i)/sigma_i
      ELSE IF (p_int .EQ. 3) THEN ! kronecker p3
        xpsij_i(ip,ij)    =+ k_T*SQRT3/SQRT2               * SQRT(tau_i)/sigma_i
        xpsijp1_i(ip,ij)  = 0._dp; xpsijm1_i(ip,ij)  = 0._dp;
      ELSE
        xpsij_i(ip,ij)    = 0._dp; xpsijp1_i(ip,ij)  = 0._dp
        xpsijm1_i(ip,ij)  = 0._dp;
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE compute_lin_coeff

!******************************************************************************!
!!!!!!! Routine that can artificially increase or wipe modes
!******************************************************************************!
SUBROUTINE save_EM_ZF_modes
  USE fields
  USE array, ONLY : moments_e_ZF, moments_i_ZF, phi_ZF, moments_e_EM,moments_i_EM,phi_EM
  USE grid
  USE time_integration, ONLY: updatetlevel
  USE initial_par, ONLY: ACT_ON_MODES
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

END MODULE numerics
