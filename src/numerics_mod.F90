!! MODULE NUMERICS
!   The module numerics contains a set of routines that are called only once at
! the begining of a run. These routines do not need to be optimzed
MODULE numerics
    USE basic
    USE prec_const
    USE grid
    USE utility

    implicit none

    PUBLIC :: build_dnjs_table, evaluate_kernels, evaluate_EM_op
    PUBLIC :: compute_lin_coeff

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
  USE model, ONLY : sigmae2_taue_o2, sigmai2_taui_o2, KIN_E
  IMPLICIT NONE
  INTEGER    :: j_int
  REAL(dp)   :: j_dp, y_, factj

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
    IF(j_int .LE. 0) THEN
      factj = 1._dp
      kernel_e(ij,iky,ikx,iz,eo) = 0._dp
    ELSE
      factj = GAMMA(j_dp+1._dp)
      kernel_e(ij,iky,ikx,iz,eo) = y_**j_int*EXP(-y_)/factj
    ENDIF
  ENDDO
  IF (ijs_e .EQ. 1) &
  kernel_e(ijgs_e,iky,ikx,iz,eo) = 0._dp
  ENDIF
  !!!!! Ion kernels !!!!!
  DO ij = ijgs_i, ijge_i
    j_int = jarray_i(ij)
    j_dp  = REAL(j_int,dp)
    y_    =  sigmai2_taui_o2 * kparray(iky,ikx,iz,eo)**2
    IF(j_int .LT. 0) THEN
      kernel_i(ij,iky,ikx,iz,eo) = 0._dp
    ELSE
      factj = GAMMA(j_dp+1._dp)
      kernel_i(ij,iky,ikx,iz,eo) = y_**j_int*EXP(-y_)/factj
    ENDIF
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
SUBROUTINE evaluate_EM_op
  IMPLICIT NONE

  CALL evaluate_poisson_op
  CALL evaluate_ampere_op

END SUBROUTINE evaluate_EM_op
!!!!!!! Evaluate inverse polarisation operator for Poisson equation
!******************************************************************************!
SUBROUTINE evaluate_poisson_op
  USE basic
  USE array, Only : kernel_e, kernel_i, inv_poisson_op, inv_pol_ion
  USE grid
  USE model, ONLY : qe2_taue, qi2_taui, KIN_E
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
    IF (KIN_E) THEN ! Kinetic model
    ! loop over n only if the max polynomial degree
    DO ine=1,jmaxe+1 ! ine = n+1
      pol_e = pol_e  + qe2_taue*kernel_e(ine,iky,ikx,iz,0)**2 ! ... sum recursively ...
    END DO
    ELSE ! Adiabatic model
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
  USE model, ONLY : q_e, q_i, beta, sigma_e, sigma_i
  USE geometry, ONLY : hatB
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
      pol_i = 0._dp
      ! loop over n only up to the max polynomial degree
      DO ini=1,jmaxi+1
        pol_i = pol_i  + kernel_i(ini,iky,ikx,iz,0)**2 ! ... sum recursively ...
      END DO
      pol_i = q_i**2/(sigma_i**2) * pol_i
      !!!!!!!!!!!!! Electron contribution
      pol_e = 0._dp
      ! loop over n only up to the max polynomial degree
      DO ine=1,jmaxe+1 ! ine = n+1
        pol_e = pol_e  + kernel_e(ine,iky,ikx,iz,0)**2 ! ... sum recursively ...
      END DO
      pol_e = q_e**2/(sigma_e**2) * pol_e
      inv_ampere_op(iky, ikx, iz) =  1._dp/(2._dp*kperp2*hatB(iz,0)**2 + beta*(pol_i + pol_e))
    ENDIF
    END DO zloop
    END DO kyloop
    END DO kxloop
  ENDIF

END SUBROUTINE evaluate_ampere_op
!******************************************************************************!

SUBROUTINE compute_lin_coeff
  USE array
  USE model, ONLY: taue_qe, taui_qi, &
                   k_Te, k_Ti, k_Ne, k_Ni, CurvB, GradB, KIN_E,&
                   tau_e, tau_i, sigma_e, sigma_i
  USE prec_const
  USE grid,  ONLY: parray_e, parray_i, jarray_e, jarray_i, &
                   ip,ij, ips_e,ipe_e, ips_i,ipe_i, ijs_e,ije_e, ijs_i,ije_i
  IMPLICIT NONE
  INTEGER     :: p_int, j_int ! polynom. degrees
  REAL(dp)    :: p_dp, j_dp
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
          xphij_e(ip,ij)    = +k_Ne+ 2.*j_dp*k_Te
          xphijp1_e(ip,ij)  = -k_Te*(j_dp+1._dp)
          xphijm1_e(ip,ij)  = -k_Te* j_dp
        ELSE IF (p_int .EQ. 2) THEN ! kronecker p2
          xphij_e(ip,ij)    = +k_Te/SQRT2
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
        xphij_i(ip,ij)    = +k_Ni + 2._dp*j_dp*k_Ti
        xphijp1_i(ip,ij)  = -k_Ti*(j_dp+1._dp)
        xphijm1_i(ip,ij)  = -k_Ti* j_dp
      ELSE IF (p_int .EQ. 2) THEN ! kronecker p2
        xphij_i(ip,ij)    = +k_Ti/SQRT2
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
          xpsij_e  (ip,ij)  = +(k_Ne + (2._dp*j_dp+1._dp)*k_Te)* SQRT(tau_e)/sigma_e
          xpsijp1_e(ip,ij)  = - k_Te*(j_dp+1._dp)              * SQRT(tau_e)/sigma_e
          xpsijm1_e(ip,ij)  = - k_Te* j_dp                     * SQRT(tau_e)/sigma_e
        ELSE IF (p_int .EQ. 3) THEN ! kronecker p3
          xpsij_e  (ip,ij)  = + k_Te*SQRT3/SQRT2               * SQRT(tau_e)/sigma_e
          xpsijp1_e(ip,ij)  = 0._dp; xpsijm1_e(ip,ij)  = 0._dp;
        ELSE
          xpsij_e  (ip,ij)  = 0._dp; xpsijp1_e(ip,ij)  = 0._dp
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
        xpsij_i  (ip,ij)  = +(k_Ni + (2._dp*j_dp+1._dp)*k_Ti)* SQRT(tau_i)/sigma_i
        xpsijp1_i(ip,ij)  = - k_Ti*(j_dp+1._dp)              * SQRT(tau_i)/sigma_i
        xpsijm1_i(ip,ij)  = - k_Ti* j_dp                     * SQRT(tau_i)/sigma_i
      ELSE IF (p_int .EQ. 3) THEN ! kronecker p3
        xpsij_i  (ip,ij)  = + k_Ti*SQRT3/SQRT2               * SQRT(tau_i)/sigma_i
        xpsijp1_i(ip,ij)  = 0._dp; xpsijm1_i(ip,ij)  = 0._dp;
      ELSE
        xpsij_i  (ip,ij)  = 0._dp; xpsijp1_i(ip,ij)  = 0._dp
        xpsijm1_i(ip,ij)  = 0._dp;
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE compute_lin_coeff

END MODULE numerics
