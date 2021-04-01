!_____________________________________________________________________________!
!_____________________________________________________________________________!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Electrons moments RHS !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!_____________________________________________________________________________!
SUBROUTINE moments_eq_rhs_e

  USE basic
  USE time_integration
  USE array
  USE fields
  USE grid
  USE model
  USE prec_const
  USE utility, ONLY : is_nan
  USE collision
  IMPLICIT NONE

  INTEGER     :: ip2, ij2, il, p_int, j_int, p2_int, j2_int ! loops indices and polynom. degrees
  REAL(dp)    :: p_dp, j_dp, l_dp
  REAL(dp)    :: kr, kz, kperp2
  REAL(dp)    :: kernelj, kerneljp1, kerneljm1 ! Kernel functions and variable
  REAL(dp)    :: xNapj, xNapp1j, xNapm1j, xNapp2j, xNapm2j, xNapjp1, xNapjm1 ! Mom. factors depending on the pj loop
  REAL(dp)    :: xphij, xphijp1, xphijm1, xphijpar ! ESpot. factors depending on the pj loop
  REAL(dp)    :: xCapj,   xCa20,   xCa01, xCa10 ! Coll. factors depending on the pj loop
  COMPLEX(dp) :: TNapj, TNapp1j, TNapm1j, TNapp2j, TNapm2j, TNapjp1, TNapjm1, Tphi
  COMPLEX(dp) :: TColl, TColl20, TColl01, TColl10 ! terms of the rhs
  COMPLEX(dp) :: i_kz, Hyper_diff_p, Hyper_diff_j

  ! Measuring execution time
  CALL cpu_time(t0_rhs)

  ploope : DO ip = ips_e, ipe_e ! loop over Hermite degree indices
    p_int= parray_e(ip) ! Hermite polynom. degree
    p_dp = REAL(p_int,dp)  ! REAL of the Hermite degree

    ! N_e^{p+1,j} coeff
    xNapp1j = sqrtTaue_qe * SQRT(p_dp + 1)
    ! N_e^{p-1,j} coeff
    xNapm1j = sqrtTaue_qe * SQRT(p_dp)

    ! N_e^{p+2,j} coeff
    xNapp2j = taue_qe_etaB * SQRT((p_dp + 1._dp) * (p_dp + 2._dp))
    ! N_e^{p-2,j} coeff
    xNapm2j = taue_qe_etaB * SQRT(p_dp * (p_dp - 1._dp))

    jloope : DO ij = ijs_e, ije_e ! loop over Laguerre degree indices
      j_int= jarray_e(ij) ! Laguerre polynom. degree
      j_dp = REAL(j_int,dp) ! REAL of degree

      ! N_e^{p,j+1} coeff
      xNapjp1 = -taue_qe_etaB * (j_dp + 1._dp)
      ! N_e^{p,j-1} coeff
      xNapjm1 = -taue_qe_etaB * j_dp
      ! N_e^{pj} coeff
      xNapj   = taue_qe_etaB * 2._dp*(p_dp + j_dp + 1._dp)

      !! Collision operator pj terms
      xCapj = -nu_e*(p_dp + 2._dp*j_dp) !DK Lenard-Bernstein basis
      ! Dougherty part
      IF ( CO .EQ. -2) THEN
        IF     ((p_int .EQ. 2) .AND. (j_int .EQ. 0)) THEN ! kronecker pj20
          xCa20 = nu_e * 2._dp/3._dp
          xCa01 = -SQRT2 * xCa20
          xCa10 = 0._dp
        ELSEIF ((p_int .EQ. 0) .AND. (j_int .EQ. 1)) THEN ! kronecker pj01
          xCa20 = -nu_e * SQRT2 * 2._dp/3._dp
          xCa01 = -SQRT2 * xCa20
          xCa10 = 0._dp
        ELSEIF ((p_int .EQ. 1) .AND. (j_int .EQ. 0)) THEN ! kronecker pj10
          xCa20 = 0._dp
          xCa01 = 0._dp
          xCa10 = nu_e
        ELSE
          xCa20 = 0._dp; xCa01 = 0._dp; xCa10 = 0._dp
        ENDIF
      ENDIF

      !! Electrostatic potential pj terms
      IF (p_int .EQ. 0) THEN ! kronecker p0
        xphij   =  (eta_n + 2.*j_dp*eta_T - 2._dp*eta_B*(j_dp+1._dp) )
        xphijp1 = -(eta_T - eta_B)*(j_dp+1._dp)
        xphijm1 = -(eta_T - eta_B)* j_dp
        xphijpar= 0._dp
      ELSE IF (p_int .EQ. 1) THEN ! kronecker p1
        xphijpar =  qe_sigmae_sqrtTaue
        xphij = 0._dp; xphijp1 = 0._dp; xphijm1 = 0._dp
      ELSE IF (p_int .EQ. 2) THEN ! kronecker p2
        xphij   =  (eta_T/SQRT2 - SQRT2*eta_B)
        xphijp1 = 0._dp; xphijm1 = 0._dp; xphijpar= 0._dp
      ELSE
        xphij = 0._dp; xphijp1 = 0._dp; xphijm1 = 0._dp; xphijpar= 0._dp
      ENDIF

      ! Loop on kspace
      krloope : DO ikr = ikrs,ikre
        kzloope : DO ikz = ikzs,ikze
          kr     = krarray(ikr)   ! Poloidal wavevector
          kz     = kzarray(ikz)   ! Toroidal wavevector
          i_kz   = imagu * kz     ! Ddz derivative
          IF (Nkz .EQ. 1) i_kz = imagu * krarray(ikr) ! If 1D simulation we put kr as kz
          kperp2 = kr**2 + kz**2  ! perpendicular wavevector

          !! Compute moments mixing terms
          ! term propto N_e^{p,j}
          TNapj    = xNapj   * moments_e(ip,ij,ikr,ikz,updatetlevel)
          ! term propto N_e^{p+2,j}
          TNapp2j  = xNapp2j * moments_e(ip+2,ij,ikr,ikz,updatetlevel)
          ! term propto N_e^{p-2,j}
          TNapm2j  = xNapm2j * moments_e(ip-2,ij,ikr,ikz,updatetlevel)
          ! term propto N_e^{p,j+1}
          TNapjp1  = xNapjp1 * moments_e(ip,ij+1,ikr,ikz,updatetlevel)
          ! term propto N_e^{p,j-1}
          TNapjm1  = xNapjm1 * moments_e(ip,ij-1,ikr,ikz,updatetlevel)

          !! Collision
          IF (CO .EQ. 0) THEN ! Lenard Bernstein
            TColl = xCapj * moments_e(ip,ij,ikr,ikz,updatetlevel)

          ELSEIF (CO .EQ. -1) THEN ! DK Dougherty
            TColl20 = 0._dp; TColl01 = 0._dp; TColl10 = 0._dp
            IF ( (pmaxe .GE. 2) ) TColl20 = xCa20 * moments_e(3,1,ikr,ikz,updatetlevel)
            IF ( (jmaxe .GE. 1) ) TColl01 = xCa01 * moments_e(1,2,ikr,ikz,updatetlevel)
            IF ( (pmaxe .GE. 1) ) TColl10 = xCa10 * moments_e(2,1,ikr,ikz,updatetlevel)
            ! Total collisional term
            TColl =  xCapj* moments_e(ip,ij,ikr,ikz,updatetlevel)&
                   + TColl20 + TColl01 + TColl10

          ELSEIF (CO .EQ. 1) THEN ! GK Dougherty
            CALL DoughertyGK_e(ip,ij,ikr,ikz,TColl)

          ELSE ! COSOLver matrix
            TColl = TColl_e(ip,ij,ikr,ikz)
        ENDIF

          !! Electrical potential term
          IF ( p_int .LE. 2 ) THEN ! kronecker p0 p1 p2
            Tphi = phi(ikr,ikz) * (xphij*kernel_e(ij, ikr, ikz) &
                     + xphijp1*kernel_e(ij+1, ikr, ikz) &
                     + xphijm1*kernel_e(ij-1, ikr, ikz) )
          ELSE
            Tphi = 0._dp
          ENDIF

          !! Parallel kinetic hyperdiffusion (projection of d/dv^4 f on Hermite basis)
          Hyper_diff_p = 0._dp
          IF ( p_int .GE. 4 ) THEN
            Hyper_diff_p = 4._dp*SQRT(p_dp*(p_dp-1._dp)*(p_dp-2._dp)*(p_dp-3._dp)*(p_dp-4._dp))&
                          *moments_e(ip-4,ij,ikr,ikz,updatetlevel)
          ENDIF

          !! Perpendicular kinetic hyperdiffusion (projection of d/dv^4 f on Laguerre basis)
          Hyper_diff_j = 0._dp
          IF ( j_int .GE. 4 ) THEN
            DO il = 1,(ij-4)
              l_dp = real(il-1,dp)
              Hyper_diff_j = Hyper_diff_j + (j_dp-(l_dp+1_dp))*moments_e(ip,il,ikr,ikz,updatetlevel)
            ENDDO
          ENDIF

          !! Sum of all linear terms
          moments_rhs_e(ip,ij,ikr,ikz,updatetlevel) = &
              -i_kz  * (TNapj + TNapp2j + TNapm2j + TNapjp1 + TNapjm1 - Tphi)&
              - mu*kperp2**2 * moments_e(ip,ij,ikr,ikz,updatetlevel) &
              + mu_p * Hyper_diff_p + mu_j * Hyper_diff_j &
              + TColl

          !! Adding non linearity
          IF ( NON_LIN ) THEN
            moments_rhs_e(ip,ij,ikr,ikz,updatetlevel) = &
              moments_rhs_e(ip,ij,ikr,ikz,updatetlevel) + Sepj(ip,ij,ikr,ikz)
          ENDIF

        END DO kzloope
      END DO krloope

    END DO jloope
  END DO ploope

  ! Execution time end
  CALL cpu_time(t1_rhs)
  tc_rhs = tc_rhs + (t1_rhs-t0_rhs)

END SUBROUTINE moments_eq_rhs_e
!_____________________________________________________________________________!
!_____________________________________________________________________________!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Ions moments RHS !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!_____________________________________________________________________________!
SUBROUTINE moments_eq_rhs_i

  USE basic
  USE time_integration, ONLY: updatetlevel
  USE array
  USE fields
  USE grid
  USE model
  USE prec_const
  USE utility,          ONLY : is_nan
  USE collision
  IMPLICIT NONE

  INTEGER     :: ip2, ij2, il, p_int, j_int, p2_int, j2_int ! loops indices and polynom. degrees
  REAL(dp)    :: p_dp, j_dp, l_dp
  REAL(dp)    :: kr, kz, kperp2
  REAL(dp)    :: kernelj, kerneljp1, kerneljm1 ! Kernel functions and variable
  REAL(dp)    :: xNapj, xNapp1j, xNapm1j, xNapp2j, xNapm2j, xNapjp1, xNapjm1 ! Mom. factors depending on the pj loop
  REAL(dp)    :: xphij, xphijp1, xphijm1, xphijpar ! ESpot. factors depending on the pj loop
  REAL(dp)    :: xCapj,   xCa20,   xCa01, xCa10 ! Coll. factors depending on the pj loop
  COMPLEX(dp) :: TNapj, TNapp1j, TNapm1j, TNapp2j, TNapm2j, TNapjp1, TNapjm1, Tphi
  COMPLEX(dp) :: TColl, TColl20, TColl01, TColl10 ! terms of the rhs
  COMPLEX(dp) :: i_kz, Hyper_diff_p, Hyper_diff_j

  LOGICAL     :: COPY_CLOS = .false. ! To test closures
  ! LOGICAL     :: COPY_CLOS = .true. ! To test closures

  ! Measuring execution time
  CALL cpu_time(t0_rhs)

  ploopi : DO ip = ips_i, ipe_i  ! Hermite loop
    p_int= parray_i(ip)   ! Hermite degree
    p_dp = REAL(p_int,dp) ! REAL of Hermite degree

    ! N_i^{p+1,j} coeff
    xNapp1j = sqrtTaui_qi * SQRT(p_dp + 1)
    ! N_i^{p-1,j} coeff
    xNapm1j = sqrtTaui_qi * SQRT(p_dp)
    ! x N_i^{p+2,j} coeff
    xNapp2j = taui_qi_etaB * SQRT((p_dp + 1._dp) * (p_dp + 2._dp))
    ! x N_i^{p-2,j} coeff
    xNapm2j = taui_qi_etaB * SQRT(p_dp * (p_dp - 1._dp))

    jloopi : DO ij = ijs_i, ije_i  ! This loop is from 1 to jmaxi+1
      j_int= jarray_i(ij)   ! REALof Laguerre degree
      j_dp = REAL(j_int,dp) ! REALof Laguerre degree

      ! x N_i^{p,j+1} coeff
      xNapjp1 = -taui_qi_etaB * (j_dp + 1._dp)
      ! x N_i^{p,j-1} coeff
      xNapjm1 = -taui_qi_etaB * j_dp
      ! x N_i^{pj} coeff
      xNapj   = taui_qi_etaB * 2._dp*(p_dp + j_dp + 1._dp)

      !! Collision operator pj terms
      xCapj = -nu_i*(p_dp + 2._dp*j_dp) !DK Lenard-Bernstein basis
      ! Dougherty part
      IF ( CO .EQ. -2) THEN
        IF     ((p_int .EQ. 2) .AND. (j_int .EQ. 0)) THEN ! kronecker pj20
          xCa20 = nu_i * 2._dp/3._dp
          xCa01 = -SQRT2 * xCa20
          xCa10 = 0._dp
        ELSEIF ((p_int .EQ. 0) .AND. (j_int .EQ. 1)) THEN ! kronecker pj01
          xCa20 = -nu_i * SQRT2 * 2._dp/3._dp
          xCa01 = -SQRT2 * xCa20
          xCa10 = 0._dp
        ELSEIF ((p_int .EQ. 1) .AND. (j_int .EQ. 0)) THEN ! kronecker pj10
          xCa20 = 0._dp
          xCa01 = 0._dp
          xCa10 = nu_i
        ELSE
          xCa20 = 0._dp; xCa01 = 0._dp; xCa10 = 0._dp
        ENDIF
      ENDIF

      !! Electrostatic potential pj terms
      IF (p_int .EQ. 0) THEN ! krokecker p0
        xphij   =  (eta_n + 2._dp*j_dp*eta_T - 2._dp*eta_B*(j_dp+1._dp))
        xphijp1 = -(eta_T - eta_B)*(j_dp+1._dp)
        xphijm1 = -(eta_T - eta_B)* j_dp
        xphijpar =  0._dp
      ELSE IF (p_int .EQ. 1) THEN ! kronecker p1
        xphijpar =  qi_sigmai_sqrtTaui
        xphij = 0._dp; xphijp1 = 0._dp; xphijm1 = 0._dp
      ELSE IF (p_int .EQ. 2) THEN !krokecker p2
        xphij   =  (eta_T/SQRT2 - SQRT2*eta_B)
        xphijp1 = 0._dp; xphijm1 = 0._dp; xphijpar =  0._dp
      ELSE
        xphij = 0._dp; xphijp1 = 0._dp; xphijm1 = 0._dp; xphijpar =  0._dp
      ENDIF

      ! Loop on kspace
      krloopi : DO ikr = ikrs,ikre
        kzloopi : DO ikz = ikzs,ikze
          kr     = krarray(ikr)   ! Poloidal wavevector
          kz     = kzarray(ikz)   ! Toroidal wavevector
          i_kz   = imagu * kz     ! Ddz derivative
          IF (Nkz .EQ. 1) i_kz = imagu * krarray(ikr) ! If 1D simulation we put kr as kz
          kperp2 = kr**2 + kz**2  ! perpendicular wavevector

          !! Compute moments mixing terms
          ! term propto N_i^{p,j}
          TNapj   = xNapj   * moments_i(ip,ij,ikr,ikz,updatetlevel)
          ! term propto N_i^{p+2,j}
          TNapp2j = xNapp2j * moments_i(ip+2,ij,ikr,ikz,updatetlevel)
          ! term propto N_i^{p-2,j}
          TNapm2j = xNapm2j * moments_i(ip-2,ij,ikr,ikz,updatetlevel)
          ! term propto N_i^{p,j+1}
          TNapjp1 = xNapjp1 * moments_i(ip,ij+1,ikr,ikz,updatetlevel)
          ! term propto N_i^{p,j-1}
          TNapjm1 = xNapjm1 * moments_i(ip,ij-1,ikr,ikz,updatetlevel)

          !! Collision
          IF (CO .EQ. 0) THEN ! Lenard Bernstein
            TColl = xCapj * moments_i(ip,ij,ikr,ikz,updatetlevel)

          ELSEIF (CO .EQ. -1) THEN ! DK Dougherty
            TColl20 = 0._dp; TColl01 = 0._dp; TColl10 = 0._dp
            IF ( (pmaxi .GE. 2) ) TColl20 = xCa20 * moments_i(3,1,ikr,ikz,updatetlevel)
            IF ( (jmaxi .GE. 1) ) TColl01 = xCa01 * moments_i(1,2,ikr,ikz,updatetlevel)
            IF ( (pmaxi .GE. 1) ) TColl10 = xCa10 * moments_i(2,1,ikr,ikz,updatetlevel)
            ! Total collisional term
            TColl =  xCapj* moments_i(ip,ij,ikr,ikz,updatetlevel)&
                   + TColl20 + TColl01 + TColl10

          ELSEIF (CO .EQ. 1) THEN ! GK Dougherty
            CALL DoughertyGK_i(ip,ij,ikr,ikz,TColl)
          ELSE! COSOLver matrix (Sugama, Coulomb)
            TColl = TColl_i(ip,ij,ikr,ikz)
          ENDIF

          !! Electrical potential term
          IF ( p_int .LE. 2 ) THEN ! kronecker p0 p1 p2
            Tphi = phi(ikr,ikz) * (xphij*kernel_i(ij, ikr, ikz) &
                     + xphijp1*kernel_i(ij+1, ikr, ikz) &
                     + xphijm1*kernel_i(ij-1, ikr, ikz) )
          ELSE
            Tphi = 0._dp
          ENDIF

          !! Kinetic hyperdiffusion
          Hyper_diff_p = 0._dp
          IF ( p_int .GE. 4 ) THEN
            Hyper_diff_p = (2._dp*p_dp)**2 *moments_i(ip-4,ij,ikr,ikz,updatetlevel)
          ENDIF

          Hyper_diff_j = 0._dp
          IF ( j_int .GE. 2 ) THEN
            DO il = 1,(ij-2)
              l_dp = real(il-1,dp)
              Hyper_diff_j = Hyper_diff_j + (j_dp-(l_dp+1_dp))*moments_i(ip,il,ikr,ikz,updatetlevel)
            ENDDO
          ENDIF

          !! Sum of linear terms
          moments_rhs_i(ip,ij,ikr,ikz,updatetlevel) = &
              -i_kz  * (TNapj + TNapp2j + TNapm2j + TNapjp1 + TNapjm1 - Tphi)&
              - mu*kperp2**2 * moments_i(ip,ij,ikr,ikz,updatetlevel) &
              + mu_p * Hyper_diff_p + mu_j * Hyper_diff_j &
              + TColl

          !! Adding non linearity
          IF ( NON_LIN ) THEN
           moments_rhs_i(ip,ij,ikr,ikz,updatetlevel) = &
             moments_rhs_i(ip,ij,ikr,ikz,updatetlevel) + Sipj(ip,ij,ikr,ikz)
          ENDIF

        END DO kzloopi
      END DO krloopi

    END DO jloopi
  END DO ploopi

  ! Execution time end
  CALL cpu_time(t1_rhs)
  tc_rhs = tc_rhs + (t1_rhs-t0_rhs)

END SUBROUTINE moments_eq_rhs_i
