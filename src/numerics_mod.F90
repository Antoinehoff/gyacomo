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
  USE prec_const, ONLY: xp, PI
  IMPLICIT NONE
  INTEGER :: p_
  DO p_ = -2,pmax
    if (p_ < 4) THEN
      dv4_Hp_coeff(p_) = 0._xp
    ELSE
      dv4_Hp_coeff(p_) = 4_xp*SQRT(REAL((p_-3)*(p_-2)*(p_-1)*p_,xp))
    ENDIF
  ENDDO
   !we scale it w.r.t. to the max degree since
   !D_4^{v}\sim (\Delta v/2)^4 and \Delta v \sim 2pi/kvpar = pi/\sqrt{2P}
   ! dv4_Hp_coeff = dv4_Hp_coeff*(1._xp/2._xp/SQRT(REAL(pmax,xp)))**4
   dv4_Hp_coeff = dv4_Hp_coeff*(PI/2._xp/SQRT(2._xp*REAL(pmax,xp)))**4
END SUBROUTINE build_dv4Hp_table
!******************************************************************************!

!******************************************************************************!
!!!!!!! Evaluate the kernels once for all
!******************************************************************************!
SUBROUTINE evaluate_kernels
  USE basic
  USE array,   ONLY : kernel!, HF_phi_correction_operator
  USE grid,    ONLY : local_na, local_nj,ngj, local_nkx, local_nky, local_nz, ngz, jarray, kparray,&
                      nzgrid
  USE species, ONLY : sigma2_tau_o2
  USE prec_const, ONLY: xp
  IMPLICIT NONE
  INTEGER    :: j_int, ia, eo, ikx, iky, iz, ij
  REAL(xp)   :: j_xp, y_, factj

DO ia  = 1,local_na
  DO eo  = 1,nzgrid
    DO ikx = 1,local_nkx
      DO iky = 1,local_nky
        DO iz  = 1,local_nz + ngz
          DO ij = 1,local_nj + ngj
            j_int = jarray(ij)
            j_xp  = REAL(j_int,xp)
            y_    =  sigma2_tau_o2(ia) * kparray(iky,ikx,iz,eo)**2
            IF(j_int .LT. 0) THEN !ghosts values
              kernel(ia,ij,iky,ikx,iz,eo) = 0._xp
            ELSE
              factj = REAL(GAMMA(j_xp+1._xp),xp)
              kernel(ia,ij,iky,ikx,iz,eo) = y_**j_int*EXP(-y_)/factj
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  ! !! Correction term for the evaluation of the heat flux
  ! HF_phi_correction_operator(:,:,:) = &
  !        2._xp * Kernel(ia,1,:,:,:,1) &
  !       -1._xp * Kernel(ia,2,:,:,:,1)
  !
  ! DO ij = 1,local_Nj
  !   j_int = jarray(ij)
  !   j_xp  = REAL(j_int,xp)
  !   HF_phi_correction_operator(:,:,:) = HF_phi_correction_operator(:,:,:) &
  !   - Kernel(ia,ij,:,:,:,1) * (&
  !       2._xp*(j_xp+1.5_xp) * Kernel(ia,ij  ,:,:,:,1) &
  !       -     (j_xp+1.0_xp) * Kernel(ia,ij+1,:,:,:,1) &
  !       -              j_xp * Kernel(ia,ij-1,:,:,:,1))
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
  USE grid,    ONLY : local_na, local_nkx, local_nky, local_nz,&
                      kxarray, kyarray, local_nj, ngj, ngz, ieven
  USE species, ONLY : q2_tau
  USE model,   ONLY : ADIAB_E, ADIAB_I, tau_i
  USE prec_const, ONLY: xp
  IMPLICIT NONE
  REAL(xp)    :: pol_tot, operator_ion     ! (Z^2/tau (1-sum_n kernel_na^2))
  INTEGER     :: in,ikx,iky,iz,ia
  REAL(xp)    :: sumker     ! (Z_a^2/tau_a (1-sum_n kernel_na^2))

  ! This term has no staggered grid dependence. It is evalued for the
  ! even z grid since poisson uses p=0 moments and phi only.
  kxloop: DO ikx = 1,local_nkx
  kyloop: DO iky = 1,local_nky
  zloop:  DO iz  = 1,local_nz
  IF( (kxarray(ikx).EQ.0._xp) .AND. (kyarray(iky).EQ.0._xp) ) THEN
      inv_poisson_op(iky, ikx, iz) =  0._xp
      inv_pol_ion   (iky, ikx, iz) =  0._xp
  ELSE
    ! loop over n only up to the max polynomial degree
    pol_tot = 0._xp  ! total polarisation term
    a:DO ia = 1,local_na ! sum over species
    ! ia = 1
      sumker  = 0._xp  ! sum of ion polarisation term
      DO in=1,local_nj
        sumker = sumker + q2_tau(ia)*kernel(ia,in+ngj/2,iky,ikx,iz+ngz/2,ieven)**2 ! ... sum recursively ...
      END DO
      pol_tot = pol_tot + q2_tau(ia) - sumker
    ENDDO a
    operator_ion = pol_tot

    IF(ADIAB_E) & ! Adiabatic electron model
      pol_tot = pol_tot + 1._xp
      
    inv_poisson_op(iky, ikx, iz) =  1._xp/pol_tot
    inv_pol_ion   (iky, ikx, iz) =  1._xp/operator_ion

    IF(ADIAB_I) & ! adiabatic ions
      inv_poisson_op(iky, ikx, iz)  = tau_i
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
  USE prec_const,   ONLY : xp
  USE array,    ONLY : kernel, inv_ampere_op
  USE grid,     ONLY : local_na, local_nkx, local_nky, local_nz, ngz, total_nj, ngj,&
                       kparray, kxarray, kyarray, SOLVE_AMPERE, iodd
  USE model,    ONLY : beta
  USE species,  ONLY : q, sigma
  USE geometry, ONLY : hatB
  USE prec_const, ONLY: xp
  IMPLICIT NONE
  REAL(xp)    :: sum_jpol, kperp2, operator     ! (Z^2/tau (1-sum_n kernel_na^2))
  INTEGER     :: in,ikx,iky,iz,ia
  ! We do not solve Ampere if beta = 0 to spare waste of ressources
  IF(SOLVE_AMPERE) THEN
    x:DO ikx = 1,local_nkx
    y:DO iky = 1,local_nky
    z:DO iz  = 1,local_nz
    kperp2 = kparray(iky,ikx,iz+ngz/2,iodd)**2
    IF( (kxarray(ikx).EQ.0._xp) .AND. (kyarray(iky).EQ.0._xp) ) THEN
        inv_ampere_op(iky, ikx, iz) =  0._xp
    ELSE
      sum_jpol = 0._xp
      a:DO ia  = 1,local_na
        ! loop over n only up to the max polynomial degree
        j:DO in=1,total_nj
          sum_jpol = sum_jpol  + q(ia)**2/(sigma(ia)**2)*kernel(ia,in+ngj/2,iky,ikx,iz+ngz/2,iodd)**2 ! ... sum recursively ...
        END DO j
      END DO a
      operator = 2._xp*kperp2*hatB(iz+ngz/2,iodd)**2 + beta*sum_jpol
      inv_ampere_op(iky, ikx, iz) =  1._xp/operator
    ENDIF
    END DO z
    END DO y
    END DO x
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
  USE species, ONLY: k_T, k_N, tau, q, sqrt_tau_o_sigma
  USE model,   ONLY: k_cB, k_gB
  USE prec_const, ONLY: xp, SQRT2, SQRT3
  USE grid,  ONLY: parray, jarray, local_na, local_np, local_nj, ngj, ngp
  INTEGER     :: ia,ip,ij,p_int, j_int ! polynom. dagrees
  REAL(xp)    :: p_xp, j_xp

  !! linear coefficients for moment RHS !!!!!!!!!!
  DO ia = 1,local_na
    DO ip = 1,local_np
      p_int= parray(ip+ngp/2)   ! Hermite degree
      p_xp = REAL(p_int,xp) ! REAL of Hermite degree
      DO ij = 1,local_nj
        j_int= jarray(ij+ngj/2)   ! Laguerre degree
        j_xp = REAL(j_int,xp) ! REAL of Laguerre degree
        ! All Napj terms
        xnapj(ia,ip,ij) = tau(ia)/q(ia)*(k_cB*(2._xp*p_xp + 1._xp) &
                                        +k_gB*(2._xp*j_xp + 1._xp))
        ! Mirror force terms
        ynapp1j  (ia,ip,ij) = -sqrt_tau_o_sigma(ia) *      (j_xp+1._xp)*SQRT(p_xp+1._xp)
        ynapm1j  (ia,ip,ij) = -sqrt_tau_o_sigma(ia) *      (j_xp+1._xp)*SQRT(p_xp)
        ynapp1jm1(ia,ip,ij) = +sqrt_tau_o_sigma(ia) *              j_xp*SQRT(p_xp+1._xp)
        ynapm1jm1(ia,ip,ij) = +sqrt_tau_o_sigma(ia) *              j_xp*SQRT(p_xp)
        ! Trapping terms
        zNapm1j  (ia,ip,ij) = +sqrt_tau_o_sigma(ia) *(2._xp*j_xp+1._xp)*SQRT(p_xp)
        zNapm1jp1(ia,ip,ij) = -sqrt_tau_o_sigma(ia) *      (j_xp+1._xp)*SQRT(p_xp)
        zNapm1jm1(ia,ip,ij) = -sqrt_tau_o_sigma(ia) *              j_xp*SQRT(p_xp)
      ENDDO
    ENDDO
    DO ip = 1,local_np
      p_int= parray(ip+ngp/2)   ! Hermite degree
      p_xp = REAL(p_int,xp) ! REAL of Hermite degree
      ! Landau damping coefficients (ddz napj term)
      xnapp1j(ia,ip) = sqrt_tau_o_sigma(ia) * SQRT(p_xp+1._xp)
      xnapm1j(ia,ip) = sqrt_tau_o_sigma(ia) * SQRT(p_xp)
      ! Magnetic curvature term
      xnapp2j(ia,ip) = tau(ia)/q(ia) * k_cB * SQRT((p_xp+1._xp)*(p_xp + 2._xp))
      xnapm2j(ia,ip) = tau(ia)/q(ia) * k_cB * SQRT( p_xp       *(p_xp - 1._xp))
    ENDDO
    DO ij = 1,local_nj
      j_int= jarray(ij+ngj/2)   ! Laguerre degree
      j_xp = REAL(j_int,xp) ! REAL of Laguerre degree
      ! Magnetic gradient term
      xnapjp1(ia,ij) = -tau(ia)/q(ia) * k_gB * (j_xp + 1._xp)
      xnapjm1(ia,ij) = -tau(ia)/q(ia) * k_gB *  j_xp
    ENDDO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! ES linear coefficients for moment RHS !!!!!!!!!!
    DO ip = 1,local_np
      p_int= parray(ip+ngp/2)   ! Hermite degree
      DO ij = 1,local_nj
        j_int= jarray(ij+ngj/2)   ! REALof Laguerre degree
        j_xp = REAL(j_int,xp) ! REALof Laguerre degree
        !! Electrostatic potential pj terms
        IF (p_int .EQ. 0) THEN ! kronecker p0
          xphij  (ia,ip,ij)  = +k_N(ia) + 2._xp*j_xp*k_T(ia)
          xphijp1(ia,ip,ij)  = -k_T(ia)*(j_xp+1._xp)
          xphijm1(ia,ip,ij)  = -k_T(ia)* j_xp
        ELSE IF (p_int .EQ. 2) THEN ! kronecker p2
          xphij(ia,ip,ij)    = +k_T(ia)/SQRT2
          xphijp1(ia,ip,ij)  = 0._xp; xphijm1(ia,ip,ij)  = 0._xp;
        ELSE
          xphij  (ia,ip,ij)  = 0._xp; xphijp1(ia,ip,ij)  = 0._xp
          xphijm1(ia,ip,ij)  = 0._xp;
        ENDIF
      ENDDO
    ENDDO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Electromagnatic linear coefficients for moment RHS !!!!!!!!!!
    DO ip = 1,local_np
      p_int= parray(ip+ngp/2)   ! Hermite degree
      DO ij = 1,local_nj
        j_int= jarray(ij+ngj/2)   ! REALof Laguerre degree
        j_xp = REAL(j_int,xp) ! REALof Laguerre degree
        IF (p_int .EQ. 1) THEN ! kronecker p1
          xpsij  (ia,ip,ij)  = +(k_N(ia) + (2._xp*j_xp+1._xp)*k_T(ia))* sqrt_tau_o_sigma(ia)
          xpsijp1(ia,ip,ij)  = - k_T(ia)*(j_xp+1._xp)                 * sqrt_tau_o_sigma(ia)
          xpsijm1(ia,ip,ij)  = - k_T(ia)* j_xp                        * sqrt_tau_o_sigma(ia)
        ELSE IF (p_int .EQ. 3) THEN ! kronecker p3
          xpsij  (ia,ip,ij)  = + k_T(ia)*SQRT3/SQRT2                  * sqrt_tau_o_sigma(ia)
          xpsijp1(ia,ip,ij)  = 0._xp; xpsijm1(ia,ip,ij)  = 0._xp;
        ELSE
          xpsij  (ia,ip,ij)  = 0._xp; xpsijp1(ia,ip,ij)  = 0._xp
          xpsijm1(ia,ip,ij)  = 0._xp;
        ENDIF
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE compute_lin_coeff

END MODULE numerics
