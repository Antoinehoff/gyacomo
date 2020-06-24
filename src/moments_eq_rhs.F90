SUBROUTINE moments_eq_rhs

  USE basic
  USE time_integration
  USE array
  USE fields
  USE fourier_grid
  USE model

  use prec_const
  IMPLICIT NONE

  INTEGER     :: ip,ij, ikr,ikz ! loops indices
  REAL(dp)    :: ip_dp, ij_dp
  REAL(dp) :: kr, kz, kperp2
  REAL(dp) :: taue_qe_etaB, taui_qi_etaB
  REAL(dp) :: kernelj, kerneljp1, kerneljm1, b_e2, b_i2 ! Kernel functions and variable
  REAL(dp) :: factj, xb_e2, xb_i2 ! Auxiliary variables
  REAL(dp) :: xNapj, xNapp2j, xNapm2j, xNapjp1, xNapjm1 ! factors depending on the pj loop
  REAL(dp) :: xphij, xphijp1, xphijm1, xColl
  COMPLEX(dp) :: TNapj, TNapp2j, TNapm2j, TNapjp1, TNapjm1, Tphi
  COMPLEX(dp) :: TColl, TColl20, TColl01, TColl10 ! terms of the rhs

  !Precompute species dependant factors
  taue_qe_etaB  = tau_e/q_e * eta_B
  xb_e2         = sigma_e**2 * tau_e!/2.0 ! species dependant factor of the Kernel, squared
  taui_qi_etaB  = tau_i/q_i * eta_B
  xb_i2         = sigma_i**2 * tau_i!/2.0 ! species dependant factor of the Kernel, squared

  !!!!!!!!! Electrons moments RHS !!!!!!!!!
  Tphi          = 0 ! electrostatic potential term
  ploope : DO ip = ips_e, ipe_e ! This loop is from 1 to pmaxe+1
    ip_dp = REAL(ip-1.,dp) ! REAL index is one minus the loop index (0 to pmax)

    ! N_e^{p+2,j} multiplicator
    IF (ip+2 .LE. ipe_e) THEN
      xNapp2j = -taue_qe_etaB * SQRT((ip_dp + 1.) * (ip_dp + 2.))
    ELSE
      xNapp2j = 0.
    ENDIF

    ! N_e^{p-2,j} multiplicator
    IF (ip-2 .GE. ips_e) THEN
      xNapm2j = -taue_qe_etaB * SQRT(ip_dp*(ip_dp - 1.))
    ELSE
      xNapm2j = 0.
    ENDIF

    jloope : DO ij = ijs_e, ije_e ! This loop is from 1 to jmaxe+1
      ij_dp = REAL(ij-1.,dp) ! REAL index is one minus the loop index (0 to jmax)

      ! N_e^{p,j+1} multiplicator
      IF (ij+1 .LE. ije_e) THEN
        xNapjp1 = taue_qe_etaB * (ij_dp + 1.)
      ELSE
        xNapjp1 = 0.
      ENDIF

      ! N_e^{p,j-1} multiplicator
      IF (ij-1 .GE. ijs_e) THEN
        xNapjm1 = taue_qe_etaB * ij_dp
      ELSE
        xNapjm1 = 0.
      ENDIF

      ! N_e^{pj} multiplicator
      xNapj   = -taue_qe_etaB * 2.*(ip_dp + ij_dp + 1.)

      ! Collision operator (DK Lenard-Bernstein basis)
      xColl = ip_dp + 2.*ij_dp
      ! ... adding Dougherty terms
      IF ( (ip .EQ. 1) .AND. (ij .EQ. 2) ) THEN ! kronecker p0 * j1
        TColl01 = 2.0/3.0*(SQRT2*moments_e(3,1,ikr,ikz,updatetlevel)&
        - 2.0*moments_e(1,2,ikr,ikz,updatetlevel))
        TColl20 = 0.0; TColl10 = 0.0;
      ELSEIF ( (ip .EQ. 3) .AND. (ij .EQ. 1) ) THEN ! kronecker p2 * j0
        TColl20 = -SQRT2/3.0*(SQRT2*moments_e(3,1,ikr,ikz,updatetlevel)&
        - 2.0*moments_e(1,2,ikr,ikz,updatetlevel))
        TColl10 = 0.0; TColl01 = 0.0;
      ELSEIF ( (ip .EQ. 2) .AND. (ij .EQ. 1) ) THEN ! kronecker p1 * j0
        TColl10 = moments_e(2,1,ikr,ikz,updatetlevel)
        TColl20 = 0.0; TColl01 = 0.0;
      ELSE
        TColl10 = 0.0; TColl20 = 0.0; TColl01 = 0.0;
      ENDIF

      ! phi multiplicator for different kernel numbers
      IF (ip .EQ. 1) THEN !(kronecker delta_p^0)
        xphij   =  (eta_n + 2.*ij_dp*eta_T - 2.*eta_B*(ij_dp+1.) )
        xphijp1 = -(eta_T - eta_B)*(ij_dp+1.)
        xphijm1 = -(eta_T - eta_B)* ij_dp
        factj   = REAL(Factorial(ij-1),dp)
      ELSE IF (ip .EQ. 3) THEN !(kronecker delta_p^2)
        xphij   =  (eta_T/SQRT2 - SQRT2*eta_B)
        factj   = REAL(Factorial(ij-1),dp)
      ELSE
        xphij = 0.; xphijp1 = 0.; xphijm1 = 0.
        factj = 1
      ENDIF

      !write(*,*) '(ip,ij) = (', ip,',', ij,')'

      ! Loop on kspace
      krloope : DO ikr = ikrs,ikre
        kzloope : DO ikz = ikzs,ikze
          kr     = krarray(ikr)   ! Poloidal wavevector
          kz     = kzarray(ikz)   ! Toroidal wavevector
          kperp2 = kr**2 + kz**2  ! perpendicular wavevector
          b_e2   = kperp2 * xb_e2 ! Bessel argument

          !! Compute moments and mixing terms
          ! term propto N_e^{p,j}
          TNapj = moments_e(ip,ij,ikr,ikz,updatetlevel) * xNapj
          ! term propto N_e^{p+2,j}
          IF (ip+2 .LE. ipe_e) THEN
            TNapp2j = moments_e(ip+2,ij,ikr,ikz,updatetlevel) * xNapp2j
          ELSE
            TNapp2j = 0.
          ENDIF
          ! term propto N_e^{p-2,j}
          IF (ip-2 .GE. ips_e) THEN
            TNapm2j = moments_e(ip-2,ij,ikr,ikz,updatetlevel) * xNapm2j
          ELSE
            TNapm2j = 0.
          ENDIF
          ! xterm propto N_e^{p,j+1}
          IF (ij+1 .LE. ije_e) THEN
            TNapjp1 = moments_e(ip,ij+1,ikr,ikz,updatetlevel) * xNapjp1
          ELSE
            TNapjp1 = 0.
          ENDIF
          ! term propto N_e^{p,j-1}
          IF (ij-1 .GE. ijs_e) THEN
            TNapjm1 = moments_e(ip,ij-1,ikr,ikz,updatetlevel) * xNapjm1
          ELSE
            TNapjm1 = 0.
          ENDIF

          ! Collision term completed (DK Dougherty)
          TColl = -nu * (xColl * moments_e(ip,ij,ikr,ikz,updatetlevel) &
                          + TColl01 + TColl10 + TColl20)

          !! Electrical potential term
          Tphi = 0
          IF ( (ip .eq. 1) .or. (ip .eq. 3) ) THEN ! 0 otherwise (krokecker delta_p^0)
            kernelj    = b_e2**(ij-1) * exp(-b_e2)/factj
            kerneljp1  = kernelj * b_e2  /(ij_dp + 1.)
            kerneljm1  = kernelj * ij_dp / b_e2
            Tphi = (xphij*Kernelj + xphijp1*Kerneljp1 + xphijm1*Kerneljm1) * phi(ikr,ikz)
          ENDIF

          ! Sum of all precomputed terms
          moments_rhs_e(ip,ij,ikr,ikz,updatetlevel) = &
              imagu * kz * (TNapj + TNapp2j + TNapm2j + TNapjp1 + TNapjm1 + Tphi) + TColl

        END DO kzloope
      END DO krloope

    END DO jloope
  END DO ploope

  !!!!!!!!! Ions moments RHS !!!!!!!!!
  Tphi = 0 ! electrostatic potential term

  ploopi : DO ip = ips_i, ipe_i
    ip_dp = REAL(ip-1.,dp)

    ! x N_i^{p+2,j}
    IF (ip+2 .LE. ipe_i) THEN
      xNapp2j = -taui_qi_etaB * SQRT((ip_dp + 1.) * (ip_dp + 2.))
    ELSE
      xNapp2j = 0.
    ENDIF

    ! x N_i^{p-2,j}
    IF (ip-2 .GE. ips_i) THEN
      xNapm2j = -taui_qi_etaB * SQRT(ip_dp * (ip_dp - 1.))
    ELSE
      xNapm2j = 0.
    ENDIF

    jloopi : DO ij = ijs_i, ije_i
      ij_dp = REAL(ij-1.,dp)

      ! x N_i^{p,j+1}
      IF (ij+1 .LE. ije_i) THEN
        xNapjp1 = taui_qi_etaB * (ij_dp + 1.)
      ELSE
        xNapjp1 = 0.
      ENDIF

      ! x N_i^{p,j-1}
      IF (ij-1 .GE. ijs_i) THEN
        xNapjm1 = taui_qi_etaB * ij_dp
      ELSE
        xNapjm1 = 0.
      ENDIF

      ! x N_i^{pj}
      xNapj   = -taui_qi_etaB * 2.*(ip_dp + ij_dp + 1.)

      ! Collision term completed (DK Dougherty)
      TColl = -nu * (xColl * moments_i(ip,ij,ikr,ikz,updatetlevel) &
                      + TColl01 + TColl10 + TColl20)

      ! ... adding Dougherty terms
      IF ( (ip .EQ. 1) .AND. (ij .EQ. 2) ) THEN ! kronecker p0 * j1
        TColl01 = 2.0/3.0*(SQRT2*moments_i(3,1,ikr,ikz,updatetlevel)&
        - 2.0*moments_i(1,2,ikr,ikz,updatetlevel))
        TColl20 = 0.0; TColl10 = 0.0;
      ELSEIF ( (ip .EQ. 3) .AND. (ij .EQ. 1) ) THEN ! kronecker p2 * j0
        TColl20 = -SQRT2/3.0*(SQRT2*moments_i(3,1,ikr,ikz,updatetlevel)&
        - 2.0*moments_i(1,2,ikr,ikz,updatetlevel))
        TColl10 = 0.0; TColl01 = 0.0;
      ELSEIF ( (ip .EQ. 2) .AND. (ij .EQ. 1) ) THEN ! kronecker p1 * j0
        TColl10 = moments_i(2,1,ikr,ikz,updatetlevel)
        TColl20 = 0.0; TColl01 = 0.0;
      ELSE
        TColl10 = 0.0; TColl20 = 0.0; TColl01 = 0.0;
      ENDIF

      ! x phi
      IF (ip .EQ. 1) THEN !(krokecker delta_p^0)
        xphij   =  (eta_n + 2.*ij_dp*eta_T - 2.*eta_B*(ij_dp+1.) )
        xphijp1 = -(eta_T - eta_B)*(ij_dp+1.)
        xphijm1 = -(eta_T - eta_B)* ij_dp
        factj   = REAL(Factorial(ij-1),dp)
      ELSE IF (ip .EQ. 3) THEN !(krokecker delta_p^2)
        xphij   =  (eta_T/SQRT2 - SQRT2*eta_B)
        factj   = REAL(Factorial(ij-1),dp)
      ELSE
        xphij = 0.; xphijp1 = 0.; xphijm1 = 0.
        factj = 1.
      ENDIF

      ! Loop on kspace
      krloopi : DO ikr = ikrs,ikre
        kzloopi : DO ikz = ikzs,ikze
          kr     = krarray(ikr)   ! Poloidal wavevector
          kz     = kzarray(ikz)   ! Toroidal wavevector
          kperp2 = kr**2 + kz**2  ! perpendicular wavevector
          b_i2   = kperp2 * xb_i2 ! Bessel argument

          !! Compute moments and mixing terms
          ! term propto N_i^{p,j}
          TNapj = moments_i(ip,ij,ikr,ikz,updatetlevel) * xNapj
          ! term propto N_i^{p+2,j}
          IF (ip+2 .LE. ipe_i) THEN
            TNapp2j = moments_i(ip+2,ij,ikr,ikz,updatetlevel) * xNapp2j
          ELSE
            TNapp2j = 0.
          ENDIF
          ! term propto N_i^{p-2,j}
          IF (ip-2 .GE. ips_i) THEN
            TNapm2j = moments_i(ip-2,ij,ikr,ikz,updatetlevel) * xNapm2j
          ELSE
            TNapm2j = 0.
          ENDIF
          ! xterm propto N_i^{p,j+1}
          IF (ij+1 .LE. ije_i) THEN
            TNapjp1 = moments_i(ip,ij+1,ikr,ikz,updatetlevel) * xNapjp1
          ELSE
            TNapjp1 = 0.
          ENDIF
          ! term propto N_i^{p,j-1}
          IF (ij-1 .GE. ijs_i) THEN
            TNapjm1 = moments_i(ip,ij-1,ikr,ikz,updatetlevel) * xNapjm1
          ELSE
            TNapjm1 = 0.
          ENDIF

          ! Collision term completed (Dougherty)
          TColl = -nu* (xColl * moments_i(ip,ij,ikr,ikz,updatetlevel))

          !! Electrical potential term
          Tphi = 0
          IF ( (ip .eq. 1) .or. (ip .eq. 3) ) THEN ! 0 otherwise (krokecker delta_p^0, delta_p^2)
            kernelj    = b_i2**(ij-1) * exp(-b_i2)/factj
            kerneljp1  = kernelj * b_i2  /(ij_dp + 1.)
            kerneljm1  = kernelj * ij_dp / b_i2
            Tphi = (xphij*Kernelj + xphijp1*Kerneljp1 + xphijm1*Kerneljm1) * phi(ikr,ikz)
          ENDIF

          ! Sum of all precomputed terms
          moments_rhs_i(ip,ij,ikr,ikz,updatetlevel) = &
              imagu * kz * (TNapj + TNapp2j + TNapm2j + TNapjp1 + TNapjm1 + Tphi) + TColl

        END DO kzloopi
      END DO krloopi

    END DO jloopi
  END DO ploopi

END SUBROUTINE moments_eq_rhs
