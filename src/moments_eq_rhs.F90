SUBROUTINE moments_eq_rhs

  USE basic
  USE time_integration
  USE array
  USE fields
  USE grid
  USE model
  USE prec_const
  USE utility, ONLY : is_nan
  IMPLICIT NONE

  INTEGER     :: ip2, ij2 ! loops indices
  REAL(dp)    :: ip_dp, ij_dp
  REAL(dp)    :: kr, kz, kperp2
  REAL(dp)    :: taue_qe_etaB, taui_qi_etaB
  REAL(dp)    :: sqrtTaue_qe, sqrtTaui_qi, qe_sigmae_sqrtTaue, qi_sigmai_sqrtTaui
  REAL(dp)    :: kernelj, kerneljp1, kerneljm1, b_e2, b_i2 ! Kernel functions and variable
  REAL(dp)    :: factj, sigmae2_taue_o2, sigmai2_taui_o2 ! Auxiliary variables
  REAL(dp)    :: xNapj, xNapp1j, xNapm1j, xNapp2j, xNapm2j, xNapjp1, xNapjm1 ! Mom. factors depending on the pj loop
  REAL(dp)    :: xphij, xphijp1, xphijm1, xphijpar ! ESpot. factors depending on the pj loop
  REAL(dp)    :: xCapj,   xCa20,   xCa01, xCa10 ! Coll. factors depending on the pj loop
  COMPLEX(dp) :: TNapj, TNapp1j, TNapm1j, TNapp2j, TNapm2j, TNapjp1, TNapjm1, Tphi
  COMPLEX(dp) :: TColl, TColl20, TColl01, TColl10 ! terms of the rhs
  COMPLEX(dp) :: test_nan
  REAL(dp)    :: nu_e, nu_i, nu_ee, nu_ie ! Species collisional frequency

  !Precompute species dependant factors
  taue_qe_etaB    = tau_e/q_e * eta_B ! factor of the magnetic moment coupling
  taui_qi_etaB    = tau_i/q_i * eta_B
  sqrtTaue_qe     = sqrt(tau_e)/q_e   ! factor of parallel moment term
  sqrtTaui_qi     = sqrt(tau_i)/q_i
  qe_sigmae_sqrtTaue = q_e/sigma_e/SQRT(tau_e) ! factor of parallel phi term
  qi_sigmai_sqrtTaui = q_i/sigma_i/SQRT(tau_i)
  sigmae2_taue_o2 = sigma_e**2 * tau_e/2._dp ! factor of the Kernel argument
  sigmai2_taui_o2 = sigma_i**2 * tau_i/2._dp
  nu_e  = nu ! electron-ion collision frequency (where already multiplied by 0.532)
  nu_i  = nu * sigma_e * (tau_i)**(-3._dp/2._dp)/SQRT2 ! ion-ion collision frequ.
  nu_ee  = nu_e/SQRT2 ! e-e coll. frequ.
  nu_ie  = nu*sigma_e**2 ! i-e coll. frequ.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!! Electrons moments RHS !!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ploope : DO ip = ips_e, ipe_e ! This loop is from 1 to pmaxe+1
    ip_dp = REAL(ip-1,dp) ! REAL index is one minus the loop index (0 to pmaxe)

    ! N_e^{p+1,j} coeff
    xNapp1j = sqrtTaue_qe * SQRT(ip_dp + 1)
    ! N_e^{p-1,j} coeff
    xNapm1j = sqrtTaue_qe * SQRT(ip_dp)

    ! N_e^{p+2,j} coeff
    xNapp2j = taue_qe_etaB * SQRT((ip_dp + 1._dp) * (ip_dp + 2._dp))
    ! N_e^{p-2,j} coeff
    xNapm2j = taue_qe_etaB * SQRT(ip_dp * (ip_dp - 1._dp))

    factj = 1.0 ! Start of the recursive factorial

    jloope : DO ij = ijs_e, ije_e ! This loop is from 1 to jmaxe+1
      ij_dp = REAL(ij-1,dp) ! REAL index is one minus the loop index (0 to jmaxe)

      ! N_e^{p,j+1} coeff
      xNapjp1 = -taue_qe_etaB * (ij_dp + 1._dp)
      ! N_e^{p,j-1} coeff
      xNapjm1 = -taue_qe_etaB * ij_dp
      ! N_e^{pj} coeff
      xNapj   = taue_qe_etaB * 2._dp*(ip_dp + ij_dp + 1._dp)

      !! Collision operator pj terms
      xCapj = -nu_e*(ip_dp + 2._dp*ij_dp) !DK Lenard-Bernstein basis
      ! Dougherty part
      IF ( CO .EQ. -2) THEN
        IF     ((ip .EQ. 3) .AND. (ij .EQ. 1)) THEN ! kronecker pj20
          xCa20 = nu_e * 2._dp/3._dp
          xCa01 = -SQRT2 * xCa20
          xCa10 = 0._dp
        ELSEIF ((ip .EQ. 1) .AND. (ij .EQ. 2)) THEN ! kronecker pj01
          xCa20 = -nu_e * SQRT2 * 2._dp/3._dp
          xCa01 = -SQRT2 * xCa20
          xCa10 = 0._dp
        ELSEIF ((ip .EQ. 2) .AND. (ij .EQ. 1)) THEN ! kronecker pj10
          xCa20 = 0._dp
          xCa01 = 0._dp
          xCa10 = nu_e
        ELSE
          xCa20 = 0._dp; xCa01 = 0._dp; xCa10 = 0._dp
        ENDIF
      ENDIF

      !! Electrostatic potential pj terms
      IF (ip .EQ. 1) THEN ! kronecker p0
        xphij   =  (eta_n + 2.*ij_dp*eta_T - 2._dp*eta_B*(ij_dp+1._dp) )
        xphijp1 = -(eta_T - eta_B)*(ij_dp+1._dp)
        xphijm1 = -(eta_T - eta_B)* ij_dp
        xphijpar= 0._dp
      ELSE IF (ip .EQ. 2) THEN ! kronecker p1
        xphijpar =  qe_sigmae_sqrtTaue
        xphij = 0._dp; xphijp1 = 0._dp; xphijm1 = 0._dp
      ELSE IF (ip .EQ. 3) THEN ! kronecker p2
        xphij   =  (eta_T/SQRT2 - SQRT2*eta_B)
        xphijp1 = 0._dp; xphijm1 = 0._dp; xphijpar= 0._dp
      ELSE
        xphij = 0._dp; xphijp1 = 0._dp; xphijm1 = 0._dp; xphijpar= 0._dp
      ENDIF

      ! Recursive factorial
      IF (ij_dp .GT. 0) THEN
        factj = factj * ij_dp
      ELSE
        factj = 1._dp
      ENDIF

      ! Loop on kspace
      krloope : DO ikr = ikrs,ikre
        kzloope : DO ikz = ikzs,ikze
          kr     = krarray(ikr)   ! Poloidal wavevector
          kz     = kzarray(ikz)   ! Toroidal wavevector
          kperp2 = kr**2 + kz**2  ! perpendicular wavevector

          IF ( DK ) THEN ! Drift kinetic model
            b_e2   = 0._dp
          ELSE
            b_e2   = kperp2 * sigmae2_taue_o2 ! Bessel argument
          ENDIF

          !! Compute moments and mixing terms
          ! term propto N_e^{p,j}
          TNapj = xNapj * moments_e(ip,ij,ikr,ikz,updatetlevel)
          ! term propto N_e^{p+1,j}
          IF (ip+1 .LE. pmaxe+1) THEN ! OoB check
            TNapp1j = xNapp1j * moments_e(ip+1,ij,ikr,ikz,updatetlevel)
          ELSE
            TNapp1j = 0._dp
          ENDIF
          ! term propto N_e^{p-1,j}
          IF (ip-1 .GE. 1) THEN ! OoB check
            TNapm1j = xNapm1j * moments_e(ip-1,ij,ikr,ikz,updatetlevel)
          ELSE
            TNapm1j = 0._dp
          ENDIF
          ! term propto N_e^{p+2,j}
          IF (ip+2 .LE. pmaxe+1) THEN ! OoB check
            TNapp2j = xNapp2j * moments_e(ip+2,ij,ikr,ikz,updatetlevel)
          ELSE
            TNapp2j = 0._dp
          ENDIF
          ! term propto N_e^{p-2,j}
          IF (ip-2 .GE. 1) THEN ! OoB check
            TNapm2j = xNapm2j * moments_e(ip-2,ij,ikr,ikz,updatetlevel)
          ELSE
            TNapm2j = 0._dp
          ENDIF
          ! xterm propto N_e^{p,j+1}
          IF (ij+1 .LE. jmaxe+1) THEN ! OoB check
            TNapjp1 = xNapjp1 * moments_e(ip,ij+1,ikr,ikz,updatetlevel)
          ELSE
            TNapjp1 = 0._dp
          ENDIF
          ! term propto N_e^{p,j-1}
          IF (ij-1 .GE. 1) THEN ! OoB check
            TNapjm1 = xNapjm1 * moments_e(ip,ij-1,ikr,ikz,updatetlevel)
          ELSE
            TNapjm1 = 0._dp
          ENDIF

          !! Collision
          IF (CO .EQ. -2) THEN ! Dougherty Collision terms
            IF ( (pmaxe .GE. 2) ) THEN ! OoB check
              TColl20 = xCa20 * moments_e(3,1,ikr,ikz,updatetlevel)
            ELSE
              TColl20 = 0._dp
            ENDIF
            IF ( (jmaxe .GE. 1) ) THEN ! OoB check
              TColl01 = xCa01 * moments_e(1,2,ikr,ikz,updatetlevel)
            ELSE
              TColl01 = 0._dp
            ENDIF
            IF ( (pmaxe .GE. 1) ) THEN ! OoB check
              TColl10 = xCa10 * moments_e(2,1,ikr,ikz,updatetlevel)
            ELSE
              TColl10 = 0._dp
            ENDIF

            ! Total collisional term
            TColl =  xCapj* moments_e(ip,ij,ikr,ikz,updatetlevel)&
                   + TColl20 + TColl01 + TColl10

          ELSEIF (CO .EQ. -1) THEN ! Full Coulomb for electrons (COSOlver matrix)

            TColl = 0._dp ! Initialization

            ploopee: DO ip2 = 1,pmaxe+1 ! sum the electron-self and electron-ion test terms
              jloopee: DO ij2 = 1,jmaxe+1
                TColl = TColl + moments_e(ip2,ij2,ikr,ikz,updatetlevel) &
                   *( nu_e  * CeipjT(bare(ip-1,ij-1), bare(ip2-1,ij2-1)) &
                     +nu_ee * Ceepj (bare(ip-1,ij-1), bare(ip2-1,ij2-1)))
              ENDDO jloopee
            ENDDO ploopee

            ploopei: DO ip2 = 1,pmaxi+1 ! sum the electron-ion field terms
              jloopei: DO ij2 = 1,jmaxi+1
                TColl = TColl + moments_i(ip2,ij2,ikr,ikz,updatetlevel) &
                  *(nu_e * CeipjF(bare(ip-1,ij-1), bari(ip2-1,ij2-1)))
              END DO jloopei
            ENDDO ploopei

          ELSEIF (CO .EQ. 0) THEN ! Lenard Bernstein
            TColl = xCapj * moments_e(ip,ij,ikr,ikz,updatetlevel)
          ENDIF

          !! Electrical potential term
          IF ( (ip .LE. 3) ) THEN ! kronecker p0 p1 p2
            kernelj    = b_e2**(ij-1) * exp(-b_e2)/factj
            kerneljp1  = kernelj * b_e2  /(ij_dp + 1._dp)
            IF ( b_e2 .NE. 0 ) THEN
              kerneljm1  = kernelj * ij_dp / b_e2
            ELSE
              kerneljm1 = 0._dp
            ENDIF
            Tphi = (xphij*kernelj + xphijp1*kerneljp1 + xphijm1*kerneljm1) * phi(ikr,ikz)
          ELSE
            Tphi = 0._dp
          ENDIF

          ! Sum of all linear terms
          moments_rhs_e(ip,ij,ikr,ikz,updatetlevel) = &
              -imagu * kz * (TNapj + TNapp2j + TNapm2j + TNapjp1 + TNapjm1 - Tphi)&
              -imagu * kpar*(TNapp1j + TNapm1j + xphijpar*kernelj*phi(ikr,ikz)) &
              + TColl

          ! Adding non linearity and Hyperdiffusivity
          IF ( NON_LIN ) THEN
            moments_rhs_e(ip,ij,ikr,ikz,updatetlevel) = &
              moments_rhs_e(ip,ij,ikr,ikz,updatetlevel) - Sepj(ip,ij,ikr,ikz) &
              - mu*kperp2**2 * moments_rhs_e(ip,ij,ikr,ikz,updatetlevel)
          ENDIF

        END DO kzloope
      END DO krloope

    END DO jloope
  END DO ploope

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!! Ions moments RHS !!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ploopi : DO ip = ips_i, ipe_i  ! This loop is from 1 to pmaxi+1
    ip_dp = REAL(ip-1,dp) ! REAL index is one minus the loop index (0 to pmaxi)

    ! N_i^{p+1,j} coeff
    xNapp1j = sqrtTaui_qi * SQRT(ip_dp + 1)
    ! N_i^{p-1,j} coeff
    xNapm1j = sqrtTaui_qi * SQRT(ip_dp)

    ! x N_i^{p+2,j} coeff
    xNapp2j = taui_qi_etaB * SQRT((ip_dp + 1._dp) * (ip_dp + 2._dp))
    ! x N_i^{p-2,j} coeff
    xNapm2j = taui_qi_etaB * SQRT(ip_dp * (ip_dp - 1._dp))

    factj = 1._dp ! Start of the recursive factorial

    jloopi : DO ij = ijs_i, ije_i  ! This loop is from 1 to jmaxi+1
      ij_dp = REAL(ij-1,dp) ! REAL index is one minus the loop index (0 to jmaxi)

      ! x N_i^{p,j+1} coeff
      xNapjp1 = -taui_qi_etaB * (ij_dp + 1._dp)
      ! x N_i^{p,j-1} coeff
      xNapjm1 = -taui_qi_etaB * ij_dp
      ! x N_i^{pj} coeff
      xNapj   = taui_qi_etaB * 2._dp*(ip_dp + ij_dp + 1._dp)

      !! Collision operator pj terms
      xCapj = -nu_i*(ip_dp + 2._dp*ij_dp) !DK Lenard-Bernstein basis
      ! Dougherty part
      IF ( CO .EQ. -2) THEN
        IF     ((ip .EQ. 3) .AND. (ij .EQ. 1)) THEN ! kronecker pj20
          xCa20 = nu_i * 2._dp/3._dp
          xCa01 = -SQRT2 * xCa20
          xCa10 = 0._dp
        ELSEIF ((ip .EQ. 1) .AND. (ij .EQ. 2)) THEN ! kronecker pj01
          xCa20 = -nu_i * SQRT2 * 2._dp/3._dp
          xCa01 = -SQRT2 * xCa20
          xCa10 = 0._dp
        ELSEIF ((ip .EQ. 2) .AND. (ij .EQ. 1)) THEN
          xCa20 = 0._dp
          xCa01 = 0._dp
          xCa10 = nu_i
        ELSE
          xCa20 = 0._dp; xCa01 = 0._dp; xCa10 = 0._dp
        ENDIF
      ENDIF

      !! Electrostatic potential pj terms
      IF (ip .EQ. 1) THEN ! krokecker p0
        xphij   =  (eta_n + 2._dp*ij_dp*eta_T - 2._dp*eta_B*(ij_dp+1._dp))
        xphijp1 = -(eta_T - eta_B)*(ij_dp+1._dp)
        xphijm1 = -(eta_T - eta_B)* ij_dp
        xphijpar =  0._dp
      ELSE IF (ip .EQ. 2) THEN ! kronecker p1
        xphijpar =  qi_sigmai_sqrtTaui
        xphij = 0._dp; xphijp1 = 0._dp; xphijm1 = 0._dp
      ELSE IF (ip .EQ. 3) THEN !krokecker p2
        xphij   =  (eta_T/SQRT2 - SQRT2*eta_B)
        xphijp1 = 0._dp; xphijm1 = 0._dp; xphijpar =  0._dp
      ELSE
        xphij = 0._dp; xphijp1 = 0._dp; xphijm1 = 0._dp; xphijpar =  0._dp
      ENDIF

      ! Recursive factorial
      IF (ij_dp .GT. 0) THEN
        factj = factj * ij_dp
      ELSE
        factj = 1._dp
      ENDIF

      ! Loop on kspace
      krloopi : DO ikr = ikrs,ikre
        kzloopi : DO ikz = ikzs,ikze
          kr     = krarray(ikr)   ! Poloidal wavevector
          kz     = kzarray(ikz)   ! Toroidal wavevector
          kperp2 = kr**2 + kz**2  ! perpendicular wavevector

          IF ( DK ) THEN ! Drift kinetic model
            b_i2   = 0._dp
          ELSE
            b_i2   = kperp2 * sigmai2_taui_o2 ! Bessel argument
          ENDIF

          !! Compute moments and mixing terms
          ! term propto N_i^{p,j}
          TNapj = xNapj * moments_i(ip,ij,ikr,ikz,updatetlevel)
          ! term propto N_i^{p+1,j}
          IF (ip+1 .LE. pmaxi+1) THEN ! OoB check
            TNapp1j = xNapp1j * moments_i(ip+1,ij,ikr,ikz,updatetlevel)
          ELSE
            TNapp1j = 0._dp
          ENDIF
          ! term propto N_i^{p-1,j}
          IF (ip-1 .GE. 1) THEN ! OoB check
            TNapm1j = xNapm1j * moments_i(ip-1,ij,ikr,ikz,updatetlevel)
          ELSE
            TNapm1j = 0._dp
          ENDIF
          ! term propto N_i^{p+2,j}
          IF (ip+2 .LE. pmaxi+1) THEN ! OoB check
            TNapp2j = xNapp2j * moments_i(ip+2,ij,ikr,ikz,updatetlevel)
          ELSE
            TNapp2j = 0._dp
          ENDIF
          ! term propto N_i^{p-2,j}
          IF (ip-2 .GE. 1) THEN ! OoB check
            TNapm2j = xNapm2j * moments_i(ip-2,ij,ikr,ikz,updatetlevel)
          ELSE
            TNapm2j = 0._dp
          ENDIF
          ! xterm propto N_i^{p,j+1}
          IF (ij+1 .LE. jmaxi+1) THEN ! OoB check
            TNapjp1 = xNapjp1 * moments_i(ip,ij+1,ikr,ikz,updatetlevel)
          ELSE
            TNapjp1 = 0._dp
          ENDIF
          ! term propto N_i^{p,j-1}
          IF (ij-1 .GE. 1) THEN ! OoB check
            TNapjm1 = xNapjm1 * moments_i(ip,ij-1,ikr,ikz,updatetlevel)
          ELSE
            TNapjm1 = 0._dp
          ENDIF

          !! Collision
          IF (CO .EQ. -2) THEN  ! Dougherty Collision terms
            IF ( (pmaxi .GE. 2) ) THEN ! OoB check
              TColl20 = xCa20 * moments_i(3,1,ikr,ikz,updatetlevel)
            ELSE
              TColl20 = 0._dp
            ENDIF
            IF ( (jmaxi .GE. 1) ) THEN ! OoB check
              TColl01 = xCa01 * moments_i(1,2,ikr,ikz,updatetlevel)
            ELSE
              TColl01 = 0._dp
            ENDIF
            IF ( (pmaxi .GE. 1) ) THEN ! OoB check
              TColl10 = xCa10 * moments_i(2,1,ikr,ikz,updatetlevel)
            ELSE
              TColl10 = 0._dp
            ENDIF

            ! Total collisional term
            TColl =  xCapj* moments_i(ip,ij,ikr,ikz,updatetlevel)&
                   + TColl20 + TColl01 + TColl10

          ELSEIF (CO .EQ. -1) THEN !!! Full Coulomb for ions (COSOlver matrix) !!!

            TColl = 0._dp ! Initialization

            ploopii: DO ip2 = 1,pmaxi+1 ! sum the ion-self and ion-electron test terms
              jloopii: DO ij2 = 1,jmaxi+1
                TColl = TColl + moments_i(ip2,ij2,ikr,ikz,updatetlevel) &
                    *( nu_ie * CiepjT(bari(ip-1,ij-1), bari(ip2-1,ij2-1)) &
                      +nu_i  * Ciipj (bari(ip-1,ij-1), bari(ip2-1,ij2-1)))
              ENDDO jloopii
            ENDDO ploopii

            ploopie: DO ip2 = 1,pmaxe+1 ! sum the ion-electron field terms
              jloopie: DO ij2 = 1,jmaxe+1
                TColl = TColl + moments_e(ip2,ij2,ikr,ikz,updatetlevel) &
                  *(nu_ie * CiepjF(bari(ip-1,ij-1), bare(ip2-1,ij2-1)))
              ENDDO jloopie
            ENDDO ploopie

          ELSEIF (CO .EQ. 0) THEN! Lenhard Bernstein
            TColl = xCapj * moments_i(ip,ij,ikr,ikz,updatetlevel)
          ENDIF

          !! Electrical potential term
          IF ( (ip .LE. 3) ) THEN ! kronecker p0 p1 p2
            kernelj    = b_i2**(ij-1) * exp(-b_i2)/factj
            kerneljp1  = kernelj * b_i2  /(ij_dp + 1._dp)
            IF ( b_i2 .NE. 0 ) THEN
              kerneljm1  = kernelj * ij_dp / b_i2
            ELSE
              kerneljm1 = 0._dp
            ENDIF
            Tphi = (xphij*kernelj + xphijp1*kerneljp1 + xphijm1*kerneljm1) * phi(ikr,ikz)
          ELSE
            Tphi = 0._dp
          ENDIF

          ! Sum of linear terms
          moments_rhs_i(ip,ij,ikr,ikz,updatetlevel) = &
              -imagu * kz * (TNapj + TNapp2j + TNapm2j + TNapjp1 + TNapjm1 - Tphi)&
              -imagu * kpar*(TNapp1j + TNapm1j + xphijpar*kernelj*phi(ikr,ikz)) &
               + TColl

          ! Adding non linearity and Hyperdiffusivity
          IF ( NON_LIN ) THEN
           moments_rhs_i(ip,ij,ikr,ikz,updatetlevel) = &
             moments_rhs_i(ip,ij,ikr,ikz,updatetlevel) - Sipj(ip,ij,ikr,ikz)&
             - mu*kperp2**2 * moments_rhs_i(ip,ij,ikr,ikz,updatetlevel)
          ENDIF
        END DO kzloopi
      END DO krloopi

    END DO jloopi
  END DO ploopi

END SUBROUTINE moments_eq_rhs
