!! MODULE NUMERICS
!   The module numerics contains a set of routines that are called only once at
! the beginng of a run. These routines do not need to be optimzed
MODULE numerics
  USE prec_const, ONLY: xp
  IMPLICIT NONE
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
  IF(pmax .GT. 0) &
   dv4_Hp_coeff = dv4_Hp_coeff*(PI/2._xp/SQRT(2._xp*REAL(pmax,xp)))**4
END SUBROUTINE build_dv4Hp_table
!******************************************************************************!

!******************************************************************************!
!!!!!!! Evaluate the kernels once for all
!******************************************************************************!
SUBROUTINE evaluate_kernels
  USE basic
  USE prec_const, ONLY: xp
  USE array,   ONLY : kernel!, HF_phi_correction_operator
  USE grid,    ONLY : local_na, local_nj,ngj, local_nkx, local_nky, local_nz, ngz, jarray, kp2array,&
                      nzgrid
  USE species, ONLY : sigma2_tau_o2
  USE model,    ONLY : KN_MODEL, ORDER
#ifdef LAPACK
  USE model,    ONLY : ORDER_NUM, ORDER_DEN
#endif
  IMPLICIT NONE
  INTEGER    :: j_int, ia, eo, ikx, iky, iz, ij
  REAL(xp)   :: j_xp, y_, factj, sigma_i
  sigma_i = 1._xp ! trivial singe sigma_a = sqrt(m_a/m_i)

  SELECT CASE (KN_MODEL)
  CASE('taylor') ! developped with Leonhard Driever
    ! Kernels based on the ORDER_NUM order Taylor series of J0
    WRITE (*,*) 'Kernel approximation uses Taylor series with ', ORDER, ' powers of k'
    DO ia  = 1,local_na
      DO eo  = 1,nzgrid
        DO ikx = 1,local_nkx
          DO iky = 1,local_nky
            DO iz  = 1,local_nz + ngz
              DO ij = 1,local_nj + ngj
                y_    =  sigma2_tau_o2(ia) * kp2array(iky,ikx,iz,eo)
                j_int = jarray(ij)
                IF (j_int > ORDER .OR. j_int < 0) THEN
                  kernel(ia,ij,ikx,iky,iz,eo) = 0._xp
                ELSE
                  kernel(ia,ij,ikx,iky,iz,eo) = taylor_kernel_n(ORDER, j_int, y_)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  CASE ('pade')
#ifdef LAPACK
    ! Kernels based on the ORDER_NUM / ORDER_DEN Pade approximation of the kernels
    WRITE (*,*) 'Kernel approximation uses ', ORDER_NUM ,'/', ORDER_DEN, ' Pade approximation'
    DO ia  = 1,local_na
      DO eo  = 1,nzgrid
        DO ikx = 1,local_nkx
          DO iky = 1,local_nky
            DO iz  = 1,local_nz + ngz
              DO ij = 1,local_nj + ngj
                y_    =  sigma2_tau_o2(ia) * kp2array(iky,ikx,iz,eo)
                j_int = jarray(ij)
                IF (j_int > ORDER_NUM .OR. j_int < 0) THEN
                  kernel(ia,ij,ikx,iky,iz,eo) = 0._xp
                ELSE
                  kernel(ia,ij,ikx,iky,iz,eo) = pade_kernel_n(j_int, y_,ORDER_NUM,ORDER_DEN)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
#else
    error stop "ERROR STOP: Pade kernels cannot be used when LAPACK is not included (marconi?)"
#endif
  CASE DEFAULT
    DO ia  = 1,local_na
      DO eo  = 1,nzgrid
        DO ikx = 1,local_nkx
          DO iky = 1,local_nky
            DO iz  = 1,local_nz + ngz
              DO ij = 1,local_nj + ngj
                j_int = jarray(ij)
                j_xp  = REAL(j_int,xp)
                y_    =  sigma2_tau_o2(ia) * kp2array(iky,ikx,iz,eo)
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
    ENDDO
  END SELECT

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
  USE model,   ONLY : ADIAB_E, ADIAB_I, tau_i, q_i
  USE prec_const, ONLY: xp
  IMPLICIT NONE
  REAL(xp)    :: pol_tot, operator_ion
  INTEGER     :: in,ikx,iky,iz,ia
  REAL(xp)    :: sumker
  ! This term has no staggered grid dependence. It is evalued for the
  ! even z grid since poisson uses p=0 moments and phi only.
  kxloop: DO ikx = 1,local_nkx
  kyloop: DO iky = 1,local_nky
  zloop:  DO iz  = 1,local_nz
  IF( (kxarray(iky,ikx).EQ.0._xp) .AND. (kyarray(iky).EQ.0._xp) ) THEN
      inv_poisson_op(iky, ikx, iz) =  0._xp
      inv_pol_ion   (iky, ikx, iz) =  0._xp
  ELSE
    ! loop over n only up to the max polynomial degree
    pol_tot = 0._xp  ! total polarisation term
    a:DO ia = 1,local_na ! sum over species
    ! ia = 1
      sumker  = 0._xp  ! sum of ion polarisation term (Z_a^2/tau_a (1-sum_n kernel_na^2))
      DO in=1,local_nj
        sumker = sumker + q2_tau(ia)*kernel(ia,in+ngj/2,iky,ikx,iz+ngz/2,ieven)**2 ! ... sum recursively ...
      END DO
      pol_tot = pol_tot + q2_tau(ia) - sumker
    ENDDO a
    operator_ion = pol_tot

    IF(ADIAB_E) & ! Adiabatic electron model
      pol_tot = pol_tot + 1._xp

    IF(ADIAB_I) & ! adiabatic ions model, kernel_i = 0 and -q_i/tau_i*phi = rho_i
      pol_tot = pol_tot + q_i**2/tau_i
      
    inv_poisson_op(iky, ikx, iz) =  1._xp/pol_tot
    inv_pol_ion   (iky, ikx, iz) =  1._xp/operator_ion
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
                       kp2array, kxarray, kyarray, SOLVE_AMPERE, iodd
  USE model,    ONLY : beta, ADIAB_I
  USE species,  ONLY : q, sigma
  USE geometry, ONLY : hatB
  USE prec_const, ONLY: xp
  IMPLICIT NONE
  REAL(xp)    :: sum_jpol, kperp2, operator, q_i, sigma_i
  INTEGER     :: in,ikx,iky,iz,ia
  q_i     = 1._xp ! single charge ion
  sigma_i = 1._xp ! trivial singe sigma_a = sqrt(m_a/m_i)
  ! We do not solve Ampere if beta = 0 to spare waste of ressources
  IF(SOLVE_AMPERE) THEN
    x:DO ikx = 1,local_nkx
    y:DO iky = 1,local_nky
    z:DO iz  = 1,local_nz
    kperp2 = kp2array(iky,ikx,iz+ngz/2,iodd)
    IF( (kxarray(iky,ikx).EQ.0._xp) .AND. (kyarray(iky).EQ.0._xp) ) THEN
        inv_ampere_op(iky, ikx, iz) =  0._xp
    ELSE
      sum_jpol = 0._xp
      a:DO ia  = 1,local_na
        ! loop over n only up to the max polynomial degree
        DO in=1,total_nj
          sum_jpol = sum_jpol  + q(ia)**2/(sigma(ia)**2)*kernel(ia,in+ngj/2,iky,ikx,iz+ngz/2,iodd)**2 ! ... sum recursively ...
        END DO 
      END DO a
      IF(ADIAB_I) THEN 
        ! no ion contribution
      ENDIF
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
  USE model,   ONLY: k_cB, k_gB, k_mB, k_tB, k_ldB
  USE prec_const, ONLY: xp, SQRT2, SQRT3
  USE grid,  ONLY: parray, jarray, local_na, local_np, local_nj, ngj_o2,ngp_o2
  INTEGER     :: ia,ip,ij,p_int, j_int ! polynom. dagrees
  REAL(xp)    :: p_xp, j_xp

  !! linear coefficients for moment RHS !!!!!!!!!!
  DO ia = 1,local_na
    DO ip = 1,local_np
      p_int= parray(ip+ngp_o2)   ! Hermite degree
      p_xp = REAL(p_int,xp) ! REAL of Hermite degree
      DO ij = 1,local_nj
        j_int= jarray(ij+ngj_o2)   ! Laguerre degree
        j_xp = REAL(j_int,xp) ! REAL of Laguerre degree
        ! All Napj terms related to magn. curvature and perp. gradient
        xnapj(ia,ip,ij) = tau(ia)/q(ia)*(k_cB*(2._xp*p_xp + 1._xp) &
                                        +k_gB*(2._xp*j_xp + 1._xp))
        ! Mirror force terms
        ynapp1j  (ia,ip,ij) = -sqrt_tau_o_sigma(ia) *      (j_xp+1._xp)*SQRT(p_xp+1._xp) * k_mB
        ynapm1j  (ia,ip,ij) = -sqrt_tau_o_sigma(ia) *      (j_xp+1._xp)*SQRT(p_xp)       * k_mB
        ynapp1jm1(ia,ip,ij) = +sqrt_tau_o_sigma(ia) *              j_xp*SQRT(p_xp+1._xp) * k_mB
        ynapm1jm1(ia,ip,ij) = +sqrt_tau_o_sigma(ia) *              j_xp*SQRT(p_xp)       * k_mB
        ! Trapping terms
        zNapm1j  (ia,ip,ij) = +sqrt_tau_o_sigma(ia) *(2._xp*j_xp+1._xp)*SQRT(p_xp) * k_tB
        zNapm1jp1(ia,ip,ij) = -sqrt_tau_o_sigma(ia) *      (j_xp+1._xp)*SQRT(p_xp) * k_tB
        zNapm1jm1(ia,ip,ij) = -sqrt_tau_o_sigma(ia) *              j_xp*SQRT(p_xp) * k_tB
      ENDDO
    ENDDO
    DO ip = 1,local_np
      p_int= parray(ip+ngp_o2)   ! Hermite degree
      p_xp = REAL(p_int,xp) ! REAL of Hermite degree
      ! Landau damping coefficients (ddz napj term)
      xnapp1j(ia,ip) = sqrt_tau_o_sigma(ia) * SQRT(p_xp+1._xp) * k_ldB
      xnapm1j(ia,ip) = sqrt_tau_o_sigma(ia) * SQRT(p_xp)       * k_ldB
      ! Magnetic curvature term
      xnapp2j(ia,ip) = tau(ia)/q(ia) * SQRT((p_xp+1._xp)*(p_xp + 2._xp)) * k_cB
      xnapm2j(ia,ip) = tau(ia)/q(ia) * SQRT( p_xp       *(p_xp - 1._xp)) * k_cB
    ENDDO
    DO ij = 1,local_nj
      j_int= jarray(ij+ngj_o2)   ! Laguerre degree
      j_xp = REAL(j_int,xp) ! REAL of Laguerre degree
      ! Magnetic perp. gradient term
      xnapjp1(ia,ij) = -tau(ia)/q(ia) * (j_xp + 1._xp) * k_gB
      xnapjm1(ia,ij) = -tau(ia)/q(ia) *  j_xp          * k_gB
    ENDDO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! ES linear coefficients for moment RHS !!!!!!!!!!
    DO ip = 1,local_np
      p_int= parray(ip+ngp_o2)   ! Hermite degree
      DO ij = 1,local_nj
        j_int= jarray(ij+ngj_o2)   ! REALof Laguerre degree
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
      p_int= parray(ip+ngp_o2)   ! Hermite degree
      DO ij = 1,local_nj
        j_int= jarray(ij+ngj_o2)   ! REALof Laguerre degree
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


!******************************************************************************!
!!!!!!! Auxilliary kernel functions/subroutines (developped with Leonhard Driever)
!******************************************************************************!
REAL(xp) FUNCTION taylor_kernel_n(order, n, y)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: order
  INTEGER, INTENT(IN) :: n
  REAL(xp), INTENT(IN) :: y
  REAL(xp) :: sum_variable
  INTEGER  ::  m
  REAL(xp) :: m_dp, n_dp

  n_dp = REAL(n, xp)
  sum_variable = 0._xp

  DO m = n, order
      m_dp = REAL(m, xp)
      sum_variable = sum_variable + (-1._xp)**(m - n) * y**m / (GAMMA(n_dp + 1._xp) * GAMMA(m_dp - n_dp + 1._xp)) ! Denominator of m C n
    END DO
  taylor_kernel_n = sum_variable
END FUNCTION taylor_kernel_n

#ifdef LAPACK
REAL(xp) FUNCTION pade_kernel_n(n, y, N_NUM, N_DEN)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, N_NUM, N_DEN
  REAL(xp), INTENT(IN) :: y
  REAL(xp) :: pade_numerator_coeffs(N_NUM + 1), pade_denominator_coeffs(N_NUM + 1)
  REAL(xp) :: numerator_sum
  REAL(xp) :: denominator_sum
  INTEGER  :: m
  ! If N_NUM == 0, then the approximation should be the same as the Taylor approx. of N_NUM:
  IF (N_NUM == 0) THEN
      pade_kernel_n = taylor_kernel_n(N_NUM, n, y)
  ELSE
    CALL find_pade_coefficients(pade_numerator_coeffs, pade_denominator_coeffs, n, N_NUM, N_DEN)
    numerator_sum = 0
    denominator_sum = 0
    DO m = 0, N_NUM
      numerator_sum = numerator_sum + pade_numerator_coeffs(m + 1) * y ** m
    END DO
    DO m = 0, N_NUM
      denominator_sum = denominator_sum + pade_denominator_coeffs(m + 1) * y ** m
    END DO
    pade_kernel_n = numerator_sum / denominator_sum
  END IF
END FUNCTION pade_kernel_n


SUBROUTINE find_pade_coefficients(pade_num_coeffs, pade_denom_coeffs, n, N_NUM, N_DEN)
  IMPLICIT NONE
#ifdef SINGLE_PRECISION
  EXTERNAL :: SGESV ! Use DGESV rather than SGESV for double precision
#else
  EXTERNAL :: DGESV ! Use DGESV rather than SGESV for double precision
#endif
  INTEGER, INTENT (IN) ::  n, N_NUM, N_DEN ! index of the considered Kernel
	REAL(xp), INTENT(OUT) :: pade_num_coeffs(N_NUM + 1), pade_denom_coeffs(N_NUM + 1) ! OUT rather than INOUT to make sure no information is retained from previous Kernel computations
	
	REAL(xp) :: taylor_kernel_coeffs(N_NUM + N_NUM + 1), denom_matrix(N_NUM, N_NUM)
	INTEGER  :: m, j
	REAL(xp) :: m_dp, n_dp
	INTEGER  :: return_code  ! for DGESV
	REAL(xp) :: pivot(N_NUM) ! for DGESV

	n_dp = REAL(n, xp)
	
	! First find the kernel Taylor expansion coefficients
	DO m = 0, N_NUM + N_NUM ! m here counts the order of the derivatives
	   m_dp = REAL(m, xp)

	   IF (m < n) THEN
	      taylor_kernel_coeffs(m + 1) = 0
	   ELSE   
	      taylor_kernel_coeffs(m + 1) = (-1)**(n + m) / (GAMMA(n_dp + 1._xp) * GAMMA(m_dp - n_dp + 1._xp))
	   END IF
	END DO

	! Next construct the denominator solving matrix
	DO m = 1, N_NUM
	   DO j = 1, N_NUM
	      IF (N_NUM + m - j < 0) THEN
		denom_matrix(m, j) = 0
	      ELSE
	      	denom_matrix(m, j) = taylor_kernel_coeffs(N_NUM + m - j + 1)
	      END IF
	   END DO
	END DO

	! Then solve for the denominator coefficients, setting the first one to 1
	!!!! SOLVER NOT YET IMPLEMENTED!!!!
	pade_denom_coeffs(1) = 1
	pade_denom_coeffs(2:) = - taylor_kernel_coeffs(N_NUM + 2 : N_NUM + N_NUM + 1) ! First acts as RHS vector for equation, is then transformed to solution by DGESV
#ifdef SINGLE_PRECISION
  CALL SGESV(N_NUM, 1, denom_matrix, N_NUM, pivot, pade_denom_coeffs(2:), N_NUM, return_code) ! LAPACK solver for matrix equation. Note that denom_matrix is now no longer as expected
#else
  CALL DGESV(N_NUM, 1, denom_matrix, N_NUM, pivot, pade_denom_coeffs(2:), N_NUM, return_code) ! LAPACK solver for matrix equation. Note that denom_matrix is now no longer as expected
#endif


	! Print an error message in case there was a problem with the solver
	IF (return_code /= 0) THEN
	   WRITE (*,*) 'An error occurred in the solving for the Pade denominator coefficients. The error code is: ', return_code
	END IF

	! Finally compute the numerator coefficients
	DO m = 1, N_NUM + 1
	   pade_num_coeffs(m) = 0 ! As the array is not automatically filled with zeros
	   DO j = 1, m	
      	      !num_matrix(m, j) = taylor_kernel_coeffs(m - j + 1)
	      pade_num_coeffs(m) = pade_num_coeffs(m) + pade_denom_coeffs(j) * taylor_kernel_coeffs(m - j + 1)
	   END DO
	END DO

END SUBROUTINE find_pade_coefficients
#endif
!******************************************************************************!


END MODULE numerics
