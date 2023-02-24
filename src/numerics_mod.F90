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
!!!!!!! Build the fourth derivative Hermite coefficient table
!******************************************************************************!
SUBROUTINE build_dv4Hp_table
  USE array, ONLY: dv4_Hp_coeff
  USE grid, ONLY: pmaxi, pmaxe
  IMPLICIT NONE
  INTEGER :: p_, pmax_
  pmax_ = MAX(pmaxi,pmaxe)
  DO p_ = -2,pmax_
    if (p_ < 4) THEN
      dv4_Hp_coeff(p_) = 0._dp
    ELSE
      dv4_Hp_coeff(p_) = 4_dp*SQRT(REAL((p_-3)*(p_-2)*(p_-1)*p_,dp))
    ENDIF
  ENDDO
   !we scale it w.r.t. to the max degree since
   !D_4^{v}\sim (\Delta v/2)^4 and \Delta v \sim 2pi/kvpar = pi/\sqrt{2P}
   ! dv4_Hp_coeff = dv4_Hp_coeff*(1._dp/2._dp/SQRT(REAL(pmax_,dp)))**4
   dv4_Hp_coeff = dv4_Hp_coeff*(PI/2._dp/SQRT(2._dp*REAL(pmax_,dp)))**4
END SUBROUTINE build_dv4Hp_table
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
    IF(j_int .LT. 0) THEN
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
    inv_poisson_op(iky, ikx, iz) =  1._dp/(qi2_taui - pol_i + qe2_taue - pol_e)
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

  USE array, ONLY:  xnepj, &
                    ynepp1j, ynepm1j, ynepp1jm1, ynepm1jm1,&
                    zNepm1j, zNepm1jp1, zNepm1jm1,&
                    xnepp1j, xnepm1j, xnepp2j, xnepm2j,&
                    xnepjp1, xnepjm1,&
                    xphij_e, xphijp1_e, xphijm1_e,&
                    xpsij_e, xpsijp1_e, xpsijm1_e,&
                    xnipj, &
                    ynipp1j, ynipm1j, ynipp1jm1, ynipm1jm1,&
                    zNipm1j, zNipm1jp1, zNipm1jm1,&
                    xnipp1j, xnipm1j, xnipp2j, xnipm2j,&
                    xnipjp1, xnipjm1,&
                    xphij_i, xphijp1_i, xphijm1_i,&
                    xpsij_i, xpsijp1_i, xpsijm1_i
  USE model, ONLY: k_Te, k_Ti, k_Ne, k_Ni, k_cB, k_gB, KIN_E,&
                   tau_e, tau_i, sigma_e, sigma_i, q_e, q_i
  USE prec_const
  USE grid,  ONLY: parray_e, parray_i, jarray_e, jarray_i, &
                   ip,ij, ips_e,ipe_e, ips_i,ipe_i, ijs_e,ije_e, ijs_i,ije_i

  IF(KIN_E) THEN
  CALL lin_coeff(k_Te,k_Ne,k_cB,k_gB,tau_e,q_e,sigma_e,&
                 parray_e(ips_e:ipe_e),jarray_e(ijs_e:ije_e),ips_e,ipe_e,ijs_e,ije_e,&
                 xnepj,xnepp1j,xnepm1j,xnepp2j,xnepm2j,xnepjp1,xnepjm1,&
                 ynepp1j,ynepm1j,ynepp1jm1,ynepm1jm1,zNepm1j,zNepm1jp1,zNepm1jm1,&
                 xphij_e,xphijp1_e,xphijm1_e,xpsij_e,xpsijp1_e,xpsijm1_e)
  ENDIF

  CALL lin_coeff(k_Ti,k_Ni,k_cB,k_gB,tau_i,q_i,sigma_i,&
             parray_i(ips_i:ipe_i),jarray_i(ijs_i:ije_i),ips_i,ipe_i,ijs_i,ije_i,&
             xnipj,xnipp1j,xnipm1j,xnipp2j,xnipm2j,xnipjp1,xnipjm1,&
             ynipp1j,ynipm1j,ynipp1jm1,ynipm1jm1,zNipm1j,zNipm1jp1,zNipm1jm1,&
             xphij_i,xphijp1_i,xphijm1_i,xpsij_i,xpsijp1_i,xpsijm1_i)

  CONTAINS
    SUBROUTINE lin_coeff(k_Ta,k_Na,k_cB,k_gB,tau_a,q_a,sigma_a,&
                         parray_a,jarray_a,ips_a,ipe_a,ijs_a,ije_a,&
                         xnapj,xnapp1j,xnapm1j,xnapp2j,xnapm2j,xnapjp1,xnapjm1,&
                         ynapp1j,ynapm1j,ynapp1jm1,ynapm1jm1,zNapm1j,zNapm1jp1,zNapm1jm1,&
                         xphij_a,xphijp1_a,xphijm1_a,xpsij_a,xpsijp1_a,xpsijm1_a)
      IMPLICIT NONE
      ! INPUTS
      REAL(dp),                         INTENT(IN) :: k_Ta,k_Na,k_cB,k_gB,tau_a,q_a,sigma_a
      INTEGER,  DIMENSION(ips_a:ipe_a), INTENT(IN) :: parray_a
      INTEGER,  DIMENSION(ijs_a:ije_a), INTENT(IN) :: jarray_a
      INTEGER,                          INTENT(IN) :: ips_a,ipe_a,ijs_a,ije_a
      ! OUTPUTS (linear coefficients used in moment_eq_rhs_mod.F90)
      REAL(dp), DIMENSION(ips_a:ipe_a,ijs_a:ije_a), INTENT(OUT) :: xnapj
      REAL(dp), DIMENSION(ips_a:ipe_a),             INTENT(OUT) :: xnapp1j, xnapm1j,   xnapp2j,   xnapm2j
      REAL(dp), DIMENSION(ijs_a:ije_a),             INTENT(OUT) :: xnapjp1, xnapjm1
      REAL(dp), DIMENSION(ips_a:ipe_a,ijs_a:ije_a), INTENT(OUT) :: ynapp1j, ynapm1j,   ynapp1jm1, ynapm1jm1
      REAL(dp), DIMENSION(ips_a:ipe_a,ijs_a:ije_a), INTENT(OUT) :: zNapm1j, zNapm1jp1, zNapm1jm1
      REAL(dp), DIMENSION(ips_a:ipe_a,ijs_a:ije_a), INTENT(OUT) :: xphij_a, xphijp1_a, xphijm1_a
      REAL(dp), DIMENSION(ips_a:ipe_a,ijs_a:ije_a), INTENT(OUT) :: xpsij_a, xpsijp1_a, xpsijm1_a
      INTEGER     :: p_int, j_int ! polynom. dagrees
      REAL(dp)    :: p_dp, j_dp
      !! linear coefficients for moment RHS !!!!!!!!!!
      DO ip = ips_a, ipe_a
        p_int= parray_a(ip)   ! Hermite degree
        p_dp = REAL(p_int,dp) ! REAL of Hermite degree
        DO ij = ijs_a, ije_a
          j_int= jarray_a(ij)   ! Laguerre degree
          j_dp = REAL(j_int,dp) ! REAL of Laguerre degree
          ! All Napj terms
          xnapj(ip,ij) = tau_a/q_a*(k_cB*(2._dp*p_dp + 1._dp) &
                                 +k_gB*(2._dp*j_dp + 1._dp))
          ! Mirror force terms
          ynapp1j  (ip,ij) = -SQRT(tau_a)/sigma_a *      (j_dp+1._dp)*SQRT(p_dp+1._dp)
          ynapm1j  (ip,ij) = -SQRT(tau_a)/sigma_a *      (j_dp+1._dp)*SQRT(p_dp)
          ynapp1jm1(ip,ij) = +SQRT(tau_a)/sigma_a *              j_dp*SQRT(p_dp+1._dp)
          ynapm1jm1(ip,ij) = +SQRT(tau_a)/sigma_a *              j_dp*SQRT(p_dp)
          ! Trapping terms
          zNapm1j  (ip,ij) = +SQRT(tau_a)/sigma_a *(2._dp*j_dp+1._dp)*SQRT(p_dp)
          zNapm1jp1(ip,ij) = -SQRT(tau_a)/sigma_a *      (j_dp+1._dp)*SQRT(p_dp)
          zNapm1jm1(ip,ij) = -SQRT(tau_a)/sigma_a *              j_dp*SQRT(p_dp)
        ENDDO
      ENDDO
      DO ip = ips_a, ipe_a
        p_int= parray_a(ip)   ! Hermite degree
        p_dp = REAL(p_int,dp) ! REAL of Hermite degree
        ! Landau damping coefficients (ddz napj term)
        xnapp1j(ip) = SQRT(tau_a)/sigma_a * SQRT(p_dp+1._dp)
        xnapm1j(ip) = SQRT(tau_a)/sigma_a * SQRT(p_dp)
        ! Magnetic curvature term
        xnapp2j(ip) = tau_a/q_a * k_cB * SQRT((p_dp+1._dp)*(p_dp + 2._dp))
        xnapm2j(ip) = tau_a/q_a * k_cB * SQRT( p_dp       *(p_dp - 1._dp))
      ENDDO
      DO ij = ijs_a, ije_a
        j_int= jarray_a(ij)   ! Laguerre degree
        j_dp = REAL(j_int,dp) ! REAL of Laguerre degree
        ! Magnetic gradient term
        xnapjp1(ij) = -tau_a/q_a * k_gB * (j_dp + 1._dp)
        xnapjm1(ij) = -tau_a/q_a * k_gB *  j_dp
      ENDDO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! ES linear coefficients for moment RHS !!!!!!!!!!
      DO ip = ips_a, ipe_a
        p_int= parray_a(ip)   ! Hermite degree
        DO ij = ijs_a, ije_a
          j_int= jarray_a(ij)   ! REALof Laguerre degree
          j_dp = REAL(j_int,dp) ! REALof Laguerre degree
          !! Electrostatic potential pj terms
          IF (p_int .EQ. 0) THEN ! kronecker p0
            xphij_a(ip,ij)    = +k_Na + 2._dp*j_dp*k_Ta
            xphijp1_a(ip,ij)  = -k_Ta*(j_dp+1._dp)
            xphijm1_a(ip,ij)  = -k_Ta* j_dp
          ELSE IF (p_int .EQ. 2) THEN ! kronecker p2
            xphij_a(ip,ij)    = +k_Ta/SQRT2
            xphijp1_a(ip,ij)  = 0._dp; xphijm1_a(ip,ij)  = 0._dp;
          ELSE
            xphij_a(ip,ij)    = 0._dp; xphijp1_a(ip,ij)  = 0._dp
            xphijm1_a(ip,ij)  = 0._dp;
          ENDIF
        ENDDO
      ENDDO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Electromagnatic linear coefficients for moment RHS !!!!!!!!!!
      DO ip = ips_a, ipe_a
        p_int= parray_a(ip)   ! Hermite degree
        DO ij = ijs_a, ije_a
          j_int= jarray_a(ij)   ! REALof Laguerre degree
          j_dp = REAL(j_int,dp) ! REALof Laguerre degree
          IF (p_int .EQ. 1) THEN ! kronecker p1
            xpsij_a  (ip,ij)  = +(k_Na + (2._dp*j_dp+1._dp)*k_Ta)* SQRT(tau_a)/sigma_a
            xpsijp1_a(ip,ij)  = - k_Ta*(j_dp+1._dp)              * SQRT(tau_a)/sigma_a
            xpsijm1_a(ip,ij)  = - k_Ta* j_dp                     * SQRT(tau_a)/sigma_a
          ELSE IF (p_int .EQ. 3) THEN ! kronecker p3
            xpsij_a  (ip,ij)  = + k_Ta*SQRT3/SQRT2               * SQRT(tau_a)/sigma_a
            xpsijp1_a(ip,ij)  = 0._dp; xpsijm1_a(ip,ij)  = 0._dp;
          ELSE
            xpsij_a  (ip,ij)  = 0._dp; xpsijp1_a(ip,ij)  = 0._dp
            xpsijm1_a(ip,ij)  = 0._dp;
          ENDIF
        ENDDO
      ENDDO
    END SUBROUTINE lin_coeff
END SUBROUTINE compute_lin_coeff

END MODULE numerics
