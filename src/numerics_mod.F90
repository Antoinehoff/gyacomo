!! MODULE NUMERICS
!   The module numerics contains a set of routines that are called only once at
! the beginng of a run. These routines do not need to be optimzed
MODULE numerics
    implicit none

    PUBLIC :: build_dnjs_table, evaluate_kernels, evaluate_EM_op
    PUBLIC :: compute_lin_coeff, build_dv4Hp_table

CONTAINS

!******************************************************************************!
!!!!!!! Build the Laguerre-Laguerre coupling coefficient table for nonlin
!******************************************************************************!
SUBROUTINE build_dnjs_table
  USE array, ONLY : dnjs
  USE FMZM,  ONLY : TO_DP
  USE coeff, ONLY : ALL2L
  USE grid,  ONLY : jmax
  IMPLICIT NONE

  INTEGER :: in, ij, is, J
  INTEGER :: n_, j_, s_

  J = jmax

  DO in = 1,J+1 ! Nested dependent loops to make benefit from dnjs symmetrys
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
  USE array,      ONLY: dv4_Hp_coeff
  USE grid,       ONLY: pmax
  USE prec_const, ONLY: dp, PI
  IMPLICIT NONE
  INTEGER :: p_
  DO p_ = -2,pmax
    if (p_ < 4) THEN
      dv4_Hp_coeff(p_) = 0._dp
    ELSE
      dv4_Hp_coeff(p_) = 4_dp*SQRT(REAL((p_-3)*(p_-2)*(p_-1)*p_,dp))
    ENDIF
  ENDDO
   !we scale it w.r.t. to the max degree since
   !D_4^{v}\sim (\Delta v/2)^4 and \Delta v \sim 2pi/kvpar = pi/\sqrt{2P}
   ! dv4_Hp_coeff = dv4_Hp_coeff*(1._dp/2._dp/SQRT(REAL(pmax,dp)))**4
   dv4_Hp_coeff = dv4_Hp_coeff*(PI/2._dp/SQRT(2._dp*REAL(pmax,dp)))**4
END SUBROUTINE build_dv4Hp_table
!******************************************************************************!

!******************************************************************************!
!!!!!!! Evaluate the kernels once for all
!******************************************************************************!
SUBROUTINE evaluate_kernels
  USE basic
  USE array,   ONLY : kernel!, HF_phi_correction_operator
  USE grid,    ONLY : local_Na, local_Nj,Ngj, local_nkx, local_nky, local_nz, Ngz, jarray, kparray,&
                      nzgrid
  USE species, ONLY : sigma2_tau_o2
  USE prec_const, ONLY: dp
  IMPLICIT NONE
  INTEGER    :: j_int, ia, eo, ikx, iky, iz, ij
  REAL(dp)   :: j_dp, y_, factj

DO ia  = 1,local_Na
  DO eo  = 1,nzgrid
    DO ikx = 1,local_nkx
      DO iky = 1,local_nky
        DO iz  = 1,local_nz + Ngz
          DO ij = 1, local_nj + Ngj
            j_int = jarray(ij)
            j_dp  = REAL(j_int,dp)
            y_    =  sigma2_tau_o2(ia) * kparray(iky,ikx,iz,eo)**2
            IF(j_int .LT. 0) THEN
              kernel(ia,ij,iky,ikx,iz,eo) = 0._dp
            ELSE
              factj = GAMMA(j_dp+1._dp)
              kernel(ia,ij,iky,ikx,iz,eo) = y_**j_int*EXP(-y_)/factj
            ENDIF
          ENDDO
          ! IF (ijs .EQ. 1) & ! if ijs is 1, the ghost kernel has negative index
            ! kernel(ia,ijgs,iky,ikx,iz,eo) = 0._dp
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  ! !! Correction term for the evaluation of the heat flux
  ! HF_phi_correction_operator(:,:,:) = &
  !        2._dp * Kernel(ia,1,:,:,:,1) &
  !       -1._dp * Kernel(ia,2,:,:,:,1)
  !
  ! DO ij = 1, local_Nj
  !   j_int = jarray(ij)
  !   j_dp  = REAL(j_int,dp)
  !   HF_phi_correction_operator(:,:,:) = HF_phi_correction_operator(:,:,:) &
  !   - Kernel(ia,ij,:,:,:,1) * (&
  !       2._dp*(j_dp+1.5_dp) * Kernel(ia,ij  ,:,:,:,1) &
  !       -     (j_dp+1.0_dp) * Kernel(ia,ij+1,:,:,:,1) &
  !       -              j_dp * Kernel(ia,ij-1,:,:,:,1))
  ! ENDDO
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
  USE array,   ONLY : kernel, inv_poisson_op, inv_pol_ion
  USE grid,    ONLY : local_Na, local_nkx, local_nky, local_nz,&
                      kxarray, kyarray, local_nj, ngj, ngz, ieven
  USE species, ONLY : q2_tau
  USE model,   ONLY : ADIAB_E, tau_e
  USE prec_const, ONLY: dp
  IMPLICIT NONE
  REAL(dp)    :: pol_ion, pol_tot, operator, operator_ion     ! (Z^2/tau (1-sum_n kernel_na^2))
  INTEGER     :: in,ikx,iky,iz,ia

  ! This term has no staggered grid dependence. It is evalued for the
  ! even z grid since poisson uses p=0 moments and phi only.
  kxloop: DO ikx = 1,local_nkx
  kyloop: DO iky = 1,local_nky
  zloop:  DO iz  = 1,local_nz
  IF( (kxarray(ikx).EQ.0._dp) .AND. (kyarray(iky).EQ.0._dp) ) THEN
      inv_poisson_op(iky, ikx, iz) =  0._dp
  ELSE
    operator = 0._dp
    DO ia = 1,local_na ! sum over species
      pol_tot = 0._dp  ! total polarisation term
      pol_ion = 0._dp  ! sum of ion polarisation term
      ! loop over n only up to the max polynomial degree
      DO in=1,local_nj
        pol_tot = pol_tot  + q2_tau(ia)*kernel(ia,in+ngj/2,iky,ikx,iz+ngz/2,ieven)**2 ! ... sum recursively ...
        pol_ion = pol_ion  + q2_tau(ia)*kernel(ia,in+ngj/2,iky,ikx,iz+ngz/2,ieven)**2 !
      END DO
      operator = operator + q2_tau(ia) - pol_tot
    ENDDO
    operator_ion = operator
    IF(ADIAB_E) THEN ! Adiabatic model
      pol_tot = pol_tot +  1._dp/tau_e - 1._dp
    ENDIF
    inv_poisson_op(iky, ikx, iz) =  1._dp/operator
    inv_pol_ion   (iky, ikx, iz) =  1._dp/operator_ion
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
  USE prec_const,   ONLY : dp
  USE array,    ONLY : kernel, inv_ampere_op
  USE grid,     ONLY : local_Na, local_nkx, local_nky, local_nz, &
                       jmax, kparray, kxarray, kyarray, SOLVE_AMPERE
  USE model,    ONLY : beta
  USE species,  ONLY : q, sigma
  USE geometry, ONLY : hatB
  USE prec_const, ONLY: dp
  IMPLICIT NONE
  REAL(dp)    :: pol_tot, kperp2     ! (Z^2/tau (1-sum_n kernel_na^2))
  INTEGER     :: in,ikx,iky,iz,ia
  ! We do not solve Ampere if beta = 0 to spare waste of ressources
  IF(SOLVE_AMPERE) THEN
    DO ikx = 1,local_nkx
    DO iky = 1,local_nky
    DO iz  = 1,local_nz
    kperp2 = kparray(iky,ikx,iz,1)**2
    IF( (kxarray(ikx).EQ.0._dp) .AND. (kyarray(iky).EQ.0._dp) ) THEN
        inv_ampere_op(iky, ikx, iz) =  0._dp
    ELSE
      !!!!!!!!!!!!!!!!! Ion contribution
      pol_tot = 0._dp
      DO ia  = 1,local_na
        ! loop over n only up to the max polynomial degree
        DO in=1,jmax+1
          pol_tot = pol_tot  + q(ia)**2/(sigma(ia)**2)*kernel(ia,in,iky,ikx,iz,1)**2 ! ... sum recursively ...
        END DO
      END DO
      inv_ampere_op(iky, ikx, iz) =  1._dp/(2._dp*kperp2*hatB(iz,0)**2 + beta*pol_tot)
    ENDIF
    END DO
    END DO
    END DO
  ENDIF
END SUBROUTINE evaluate_ampere_op
!******************************************************************************!

SUBROUTINE compute_lin_coeff

  USE array, ONLY:  xnapj, &
                    ynapp1j, ynapm1j, ynapp1jm1, ynapm1jm1,&
                    zNapm1j, zNapm1jp1, zNapm1jm1,&
                    xnapj, xnapjp1, xnapjm1,&
                    xnapp1j, xnapm1j, xnapp2j, xnapm2j,&
                    xphij, xphijp1, xphijm1,&
                    xpsij, xpsijp1, xpsijm1
  USE species, ONLY: k_T, k_N, tau, q, sqrtTau_q, tau_q
  USE model,   ONLY: k_cB, k_gB
  USE prec_const, ONLY: dp, SQRT2, SQRT3
  USE grid,  ONLY: parray, jarray, local_na, local_np, local_nj, ngj, ngp
  INTEGER     :: ia,ip,ij,p_int, j_int ! polynom. dagrees
  REAL(dp)    :: p_dp, j_dp

  !! linear coefficients for moment RHS !!!!!!!!!!
  DO ia = 1,local_na
    DO ip = 1, local_np
      p_int= parray(ip+ngp/2)   ! Hermite degree
      p_dp = REAL(p_int,dp) ! REAL of Hermite degree
      DO ij = 1, local_nj
        j_int= jarray(ij+ngj/2)   ! Laguerre degree
        j_dp = REAL(j_int,dp) ! REAL of Laguerre degree
        ! All Napj terms
        xnapj(ia,ip,ij) = tau(ia)/q(ia)*(k_cB*(2._dp*p_dp + 1._dp) &
                               +k_gB*(2._dp*j_dp + 1._dp))
        ! Mirror force terms
        ynapp1j  (ia,ip,ij) = -sqrtTau_q(ia) *      (j_dp+1._dp)*SQRT(p_dp+1._dp)
        ynapm1j  (ia,ip,ij) = -sqrtTau_q(ia) *      (j_dp+1._dp)*SQRT(p_dp)
        ynapp1jm1(ia,ip,ij) = +sqrtTau_q(ia) *              j_dp*SQRT(p_dp+1._dp)
        ynapm1jm1(ia,ip,ij) = +sqrtTau_q(ia) *              j_dp*SQRT(p_dp)
        ! Trapping terms
        zNapm1j  (ia,ip,ij) = +sqrtTau_q(ia) *(2._dp*j_dp+1._dp)*SQRT(p_dp)
        zNapm1jp1(ia,ip,ij) = -sqrtTau_q(ia) *      (j_dp+1._dp)*SQRT(p_dp)
        zNapm1jm1(ia,ip,ij) = -sqrtTau_q(ia) *              j_dp*SQRT(p_dp)
      ENDDO
    ENDDO
    DO ip = 1, local_np
      p_int= parray(ip+ngp/2)   ! Hermite degree
      p_dp = REAL(p_int,dp) ! REAL of Hermite degree
      ! Landau damping coefficients (ddz napj term)
      xnapp1j(ia,ip) = sqrtTau_q(ia) * SQRT(p_dp+1._dp)
      xnapm1j(ia,ip) = sqrtTau_q(ia) * SQRT(p_dp)
      ! Magnetic curvature term
      xnapp2j(ia,ip) = tau_q(ia) * k_cB * SQRT((p_dp+1._dp)*(p_dp + 2._dp))
      xnapm2j(ia,ip) = tau_q(ia) * k_cB * SQRT( p_dp       *(p_dp - 1._dp))
    ENDDO
    DO ij = 1, local_nj
      j_int= jarray(ij+ngj/2)   ! Laguerre degree
      j_dp = REAL(j_int,dp) ! REAL of Laguerre degree
      ! Magnetic gradient term
      xnapjp1(ia,ij) = -tau_q(ia) * k_gB * (j_dp + 1._dp)
      xnapjm1(ia,ij) = -tau_q(ia) * k_gB *  j_dp
    ENDDO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! ES linear coefficients for moment RHS !!!!!!!!!!
    DO ip = 1, local_np
      p_int= parray(ip+ngp/2)   ! Hermite degree
      DO ij = 1, local_nj
        j_int= jarray(ij+ngj/2)   ! REALof Laguerre degree
        j_dp = REAL(j_int,dp) ! REALof Laguerre degree
        !! Electrostatic potential pj terms
        IF (p_int .EQ. 0) THEN ! kronecker p0
          xphij  (ia,ip,ij)  = +k_N(ia) + 2._dp*j_dp*k_T(ia)
          xphijp1(ia,ip,ij)  = -k_T(ia)*(j_dp+1._dp)
          xphijm1(ia,ip,ij)  = -k_T(ia)* j_dp
        ELSE IF (p_int .EQ. 2) THEN ! kronecker p2
          xphij(ia,ip,ij)    = +k_T(ia)/SQRT2
          xphijp1(ia,ip,ij)  = 0._dp; xphijm1(ia,ip,ij)  = 0._dp;
        ELSE
          xphij  (ia,ip,ij)  = 0._dp; xphijp1(ia,ip,ij)  = 0._dp
          xphijm1(ia,ip,ij)  = 0._dp;
        ENDIF
      ENDDO
    ENDDO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Electromagnatic linear coefficients for moment RHS !!!!!!!!!!
    DO ip = 1, local_np
      p_int= parray(ip+ngp/2)   ! Hermite degree
      DO ij = 1, local_nj
        j_int= jarray(ij+ngj/2)   ! REALof Laguerre degree
        j_dp = REAL(j_int,dp) ! REALof Laguerre degree
        IF (p_int .EQ. 1) THEN ! kronecker p1
          xpsij  (ia,ip,ij)  = +(k_N(ia) + (2._dp*j_dp+1._dp)*k_T(ia))* sqrtTau_q(ia)
          xpsijp1(ia,ip,ij)  = - k_T(ia)*(j_dp+1._dp)              * sqrtTau_q(ia)
          xpsijm1(ia,ip,ij)  = - k_T(ia)* j_dp                     * sqrtTau_q(ia)
        ELSE IF (p_int .EQ. 3) THEN ! kronecker p3
          xpsij  (ia,ip,ij)  = + k_T(ia)*SQRT3/SQRT2               * sqrtTau_q(ia)
          xpsijp1(ia,ip,ij)  = 0._dp; xpsijm1(ia,ip,ij)  = 0._dp;
        ELSE
          xpsij  (ia,ip,ij)  = 0._dp; xpsijp1(ia,ip,ij)  = 0._dp
          xpsijm1(ia,ip,ij)  = 0._dp;
        ENDIF
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE compute_lin_coeff

END MODULE numerics
