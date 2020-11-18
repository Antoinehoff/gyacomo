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

  INTEGER     :: ip2, ij2, p_int, j_int, p2_int, j2_int ! loops indices and polynom. degrees
  REAL(dp)    :: p_dp, j_dp
  REAL(dp)    :: kr, kz, kperp2
  REAL(dp)    :: taue_qe_etaB, taui_qi_etaB
  REAL(dp)    :: sqrtTaue_qe, sqrtTaui_qi, qe_sigmae_sqrtTaue, qi_sigmai_sqrtTaui
  REAL(dp)    :: kernelj, kerneljp1, kerneljm1, be_2, bi_2 ! Kernel functions and variable
  REAL(dp)    :: factj, sigmae2_taue_o2, sigmai2_taui_o2 ! Auxiliary variables
  REAL(dp)    :: xNapj, xNapp1j, xNapm1j, xNapp2j, xNapm2j, xNapjp1, xNapjm1 ! Mom. factors depending on the pj loop
  REAL(dp)    :: xphij, xphijp1, xphijm1, xphijpar ! ESpot. factors depending on the pj loop
  REAL(dp)    :: xCapj,   xCa20,   xCa01, xCa10 ! Coll. factors depending on the pj loop
  COMPLEX(dp) :: TNapj, TNapp1j, TNapm1j, TNapp2j, TNapm2j, TNapjp1, TNapjm1, Tphi
  COMPLEX(dp) :: TColl, TColl20, TColl01, TColl10 ! terms of the rhs
  COMPLEX(dp) :: test_nan
  REAL(dp)    :: nu_e, nu_i, nu_ee, nu_ie ! Species collisional frequency

  ! Measuring execution time
  CALL cpu_time(t0_rhs)

  !Precompute species dependant factors
  IF( q_e .NE. 0._dp ) THEN
    taue_qe_etaB    = tau_e/q_e * eta_B ! factor of the magnetic moment coupling
    sqrtTaue_qe     = sqrt(tau_e)/q_e   ! factor of parallel moment term
  ELSE
    taue_qe_etaB  = 0._dp
    sqrtTaue_qe   = 0._dp
  ENDIF

  taui_qi_etaB    = tau_i/q_i * eta_B ! factor of the magnetic moment coupling
  sqrtTaui_qi     = sqrt(tau_i)/q_i   ! factor of parallel moment term
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

    factj = 1.0 ! Start of the recursive factorial

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

      ! Recursive factorial
      IF (j_dp .GT. 0) THEN
        factj = factj * j_dp
      ELSE
        factj = 1._dp
      ENDIF

      ! Loop on kspace
      krloope : DO ikr = ikrs,ikre
        kzloope : DO ikz = ikzs,ikze
          kr     = krarray(ikr)   ! Poloidal wavevector
          kz     = kzarray(ikz)   ! Toroidal wavevector
          kperp2 = kr**2 + kz**2  ! perpendicular wavevector
          be_2   = kperp2 * sigmae2_taue_o2 ! Kernel argument

          !! Compute moments and mixing terms
          ! term propto N_e^{p,j}
          TNapj = xNapj * moments_e(ip,ij,ikr,ikz,updatetlevel)
          ! term propto N_e^{p+1,j} and kparallel
          IF ( (ip+1 .LE. pmaxe+1) .AND. (.NOT. CANCEL_ODD_P) ) THEN ! OoB check
            TNapp1j = xNapp1j * moments_e(ip+1,ij,ikr,ikz,updatetlevel)
          ELSE
            TNapp1j = 0._dp
          ENDIF
          ! term propto N_e^{p-1,j} and kparallel
          IF ( (ip-1 .GE. 1) .AND. (.NOT. CANCEL_ODD_P) ) THEN ! OoB check
            TNapm1j = xNapm1j * moments_e(ip-1,ij,ikr,ikz,updatetlevel)
          ELSE
            TNapm1j = 0._dp
          ENDIF
          ! term propto N_e^{p+2,j}
          IF (ip+(2-pskip) .LE. pmaxe+1) THEN ! OoB check
            TNapp2j = xNapp2j * moments_e(ip+(2-pskip),ij,ikr,ikz,updatetlevel)
          ELSE
            TNapp2j = 0._dp
          ENDIF
          ! term propto N_e^{p-2,j}
          IF (ip-(2-pskip) .GE. 1) THEN ! OoB check
            TNapm2j = xNapm2j * moments_e(ip-(2-pskip),ij,ikr,ikz,updatetlevel)
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
            IF ( (pmaxe .GE. 2-pskip) ) THEN ! OoB check
              TColl20 = xCa20 * moments_e(3-pskip,1,ikr,ikz,updatetlevel)
            ELSE
              TColl20 = 0._dp
            ENDIF
            IF ( (jmaxe .GE. 1) ) THEN ! OoB check
              TColl01 = xCa01 * moments_e(1,2,ikr,ikz,updatetlevel)
            ELSE
              TColl01 = 0._dp
            ENDIF
            IF ( (pmaxe .GE. 1) .AND. (.NOT. CANCEL_ODD_P) ) THEN ! OoB + odd number for Hermite degrees check
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
              p2_int = parray_e(ip2)
              jloopee: DO ij2 = 1,jmaxe+1
                j2_int = jarray_e(ij2)

                TColl = TColl + moments_e(ip2,ij2,ikr,ikz,updatetlevel) &
                   *( nu_e  * CeipjT(bare(p_int,j_int), bare(p2_int,j2_int)) &
                     +nu_ee * Ceepj (bare(p_int,j_int), bare(p2_int,j2_int)))

              ENDDO jloopee
            ENDDO ploopee

            ploopei: DO ip2 = 1,pmaxi+1 ! sum the electron-ion field terms
              p2_int = parray_i(ip2)
              jloopei: DO ij2 = 1,jmaxi+1
                j2_int = jarray_i(ij2)

                TColl = TColl + moments_i(ip2,ij2,ikr,ikz,updatetlevel) &
                  *(nu_e * CeipjF(bare(p_int,j_int), bari(p2_int,j2_int)))

              END DO jloopei
            ENDDO ploopei

          ELSEIF (CO .EQ. 0) THEN ! Lenard Bernstein
            TColl = xCapj * moments_e(ip,ij,ikr,ikz,updatetlevel)
          ENDIF

          !! Electrical potential term
          IF ( (p_int .LE. 2) ) THEN ! kronecker p0 p1 p2
            kernelj    = be_2**j_int * exp(-be_2)/factj
            kerneljp1  = kernelj * be_2  /(j_dp + 1._dp)
            IF ( be_2 .NE. 0 ) THEN
              kerneljm1  = kernelj * j_dp / be_2
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
              - mu*kperp2**2 * moments_e(ip,ij,ikr,ikz,updatetlevel) &
              + TColl

          ! Adding non linearity
          IF ( NON_LIN .OR. (A0KH .NE. 0) ) THEN
            moments_rhs_e(ip,ij,ikr,ikz,updatetlevel) = &
              moments_rhs_e(ip,ij,ikr,ikz,updatetlevel) - Sepj(ip,ij,ikr,ikz)
          ENDIF

        END DO kzloope
      END DO krloope

    END DO jloope
  END DO ploope

  !_____________________________________________________________________________!
  !_____________________________________________________________________________!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!! Ions moments RHS !!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !_____________________________________________________________________________!
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

    factj = 1._dp ! Start of the recursive factorial

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

      ! Recursive factorial
      IF (j_dp .GT. 0) THEN
        factj = factj * j_dp
      ELSE
        factj = 1._dp
      ENDIF

      ! Loop on kspace
      krloopi : DO ikr = ikrs,ikre
        kzloopi : DO ikz = ikzs,ikze
          kr     = krarray(ikr)   ! Poloidal wavevector
          kz     = kzarray(ikz)   ! Toroidal wavevector
          kperp2 = kr**2 + kz**2  ! perpendicular wavevector
          bi_2   = kperp2 * sigmai2_taui_o2 ! Kernel argument

          !! Compute moments and mixing terms
          ! term propto N_i^{p,j}
          TNapj = xNapj * moments_i(ip,ij,ikr,ikz,updatetlevel)
          ! term propto N_i^{p+1,j}
          IF ( (ip+1 .LE. pmaxi+1)  .AND. (.NOT. CANCEL_ODD_P) ) THEN ! OoB check
            TNapp1j = xNapp1j * moments_i(ip+1,ij,ikr,ikz,updatetlevel)
          ELSE
            TNapp1j = 0._dp
          ENDIF
          ! term propto N_i^{p-1,j}
          IF ( (ip-1 .GE. 1) .AND. (.NOT. CANCEL_ODD_P) ) THEN ! OoB check
            TNapm1j = xNapm1j * moments_i(ip-1,ij,ikr,ikz,updatetlevel)
          ELSE
            TNapm1j = 0._dp
          ENDIF
          ! term propto N_i^{p+2,j}
          IF (ip+(2-pskip) .LE. pmaxi+1) THEN ! OoB check
            TNapp2j = xNapp2j * moments_i(ip+(2-pskip),ij,ikr,ikz,updatetlevel)
          ELSE
            TNapp2j = 0._dp
          ENDIF
          ! term propto N_i^{p-2,j}
          IF (ip-(2-pskip) .GE. 1) THEN ! OoB check
            TNapm2j = xNapm2j * moments_i(ip-(2-pskip),ij,ikr,ikz,updatetlevel)
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
            IF ( (pmaxi .GE. 2-pskip) ) THEN ! OoB check
              TColl20 = xCa20 * moments_i(3-pskip,1,ikr,ikz,updatetlevel)
            ELSE
              TColl20 = 0._dp
            ENDIF
            IF ( (jmaxi .GE. 1) ) THEN ! OoB check
              TColl01 = xCa01 * moments_i(1,2,ikr,ikz,updatetlevel)
            ELSE
              TColl01 = 0._dp
            ENDIF
            IF ( (pmaxi .GE. 1) .AND. (.NOT. CANCEL_ODD_P) ) THEN ! OoB check
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
              p2_int = parray_i(ip2)
              jloopii: DO ij2 = 1,jmaxi+1
                j2_int = jarray_i(ij2)

                TColl = TColl + moments_i(ip2,ij2,ikr,ikz,updatetlevel) &
                    *( nu_ie * CiepjT(bari(p_int,j_int), bari(p2_int,j2_int)) &
                      +nu_i  * Ciipj (bari(p_int,j_int), bari(p2_int,j2_int)))

              ENDDO jloopii
            ENDDO ploopii

            ploopie: DO ip2 = 1,pmaxe+1 ! sum the ion-electron field terms
              p2_int = parray_e(ip2)
              jloopie: DO ij2 = 1,jmaxe+1
                j2_int = jarray_e(ij2)

                TColl = TColl + moments_e(ip2,ij2,ikr,ikz,updatetlevel) &
                  *(nu_ie * CiepjF(bari(p_int,j_int), bare(p2_int,j2_int)))

              ENDDO jloopie
            ENDDO ploopie

          ELSEIF (CO .EQ. 0) THEN! Lenhard Bernstein
            TColl = xCapj * moments_i(ip,ij,ikr,ikz,updatetlevel)
          ENDIF

          !! Electrical potential term
          IF ( (p_int .LE. 2) ) THEN ! kronecker p0 p1 p2
            kernelj    = bi_2**j_int * exp(-bi_2)/factj
            kerneljp1  = kernelj * bi_2  /(j_dp + 1._dp)
            IF ( bi_2 .NE. 0 ) THEN
              kerneljm1  = kernelj * j_dp / bi_2
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
              - mu*kperp2**2 * moments_i(ip,ij,ikr,ikz,updatetlevel) &
               + TColl

          ! Adding non linearity
          IF ( NON_LIN .OR. (A0KH .NE. 0) ) THEN
           moments_rhs_i(ip,ij,ikr,ikz,updatetlevel) = &
             moments_rhs_i(ip,ij,ikr,ikz,updatetlevel) - Sipj(ip,ij,ikr,ikz)
          ENDIF

        END DO kzloopi
      END DO krloopi

    END DO jloopi
  END DO ploopi

  ! Execution time end
  CALL cpu_time(t1_rhs)
  tc_rhs = tc_rhs + (t1_rhs-t0_rhs)

END SUBROUTINE moments_eq_rhs
