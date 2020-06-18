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
  COMPLEX(dp) :: kr, kz, kperp2 
  COMPLEX(dp) :: taue_qe_etaBi, taui_qi_etaBi
  COMPLEX(dp) :: kernelj, kerneljp1, kerneljm1, b_e2, b_i2 ! Kernel functions and variable
  COMPLEX(dp) :: factj, xb_e2, xb_i2 ! Auxiliary variables
  COMPLEX(dp) :: xNapj, xNapp2j, xNapm2j, xNapjp1, xNapjm1 ! factors depending on the pj loop
  COMPLEX(dp) :: xphij, xphijp1, xphijm1, xColl
  COMPLEX(dp) :: TNapj, TNapp2j, TNapm2j, TNapjp1, TNapjm1, Tphi, TColl ! terms of the rhs

  !!!!!!!!! Electrons moments RHS !!!!!!!!!
  Tphi    = 0 ! electrostatic potential term
  taue_qe_etaBi = tau_e/q_e * eta_B * imagu ! one-time computable factor
  xb_e2 = sigma_e**2 * tau_e/2. ! species dependant factor of the Kernel, squared

  ploope : DO ip = ips_e, ipe_e
    ip_dp = real(ip-1,dp) ! real index is one minus the loop index (0 to pmax)

    ! x N_e^{p+2,j}
    IF (ip+2 .LE. pmaxe + 1) THEN
      xNapp2j = taue_qe_etaBi * SQRT(ip_dp + 1._dp) * SQRT(ip_dp + 2._dp)
    ELSE
      xNapp2j = 0.
    ENDIF

    ! x N_e^{p-2,j}
    IF (ip-2 .GE. 1) THEN
      xNapm2j = taue_qe_etaBi * SQRT(ip_dp) * SQRT(ip_dp - 1.)
    ELSE
      xNapm2j = 0.
    ENDIF

    jloope : DO ij = ijs_e, ije_e
      ij_dp = real(ij-1,dp) ! real index is one minus the loop index (0 to jmax)

      ! x N_e^{p,j+1}
      IF (ij+1 .LE. jmaxe+1) THEN
        xNapjp1 = -taue_qe_etaBi * (ij_dp + 1.)
      ELSE
        xNapjp1 = 0.
      ENDIF

      ! x N_e^{p,j-1}
      IF (ij-1 .GE. 1) THEN
        xNapjm1 = -taue_qe_etaBi * ij_dp
      ELSE
        xNapjm1 = 0.
      ENDIF

      ! x N_e^{pj}
      xNapj   = taue_qe_etaBi * 2.*(ip_dp + ij_dp + 1.)

      ! first part of collision operator (Lenard-Bernstein)
      xColl = -nu*(ip_dp + 2.*ij_dp)

      ! x phi
      IF (ip .eq. 1) THEN !(krokecker delta_p^0)
        xphij   =  imagu * (eta_n + 2.*ij_dp*eta_T + 2.*eta_B*(ij_dp+1.) )
        xphijp1 = -imagu * (eta_B + eta_T)*(ij_dp+1.)
        xphijm1 = -imagu * (eta_B + eta_T)* ij_dp
        factj = real(Factorial(ij-1),dp)
      ELSE IF (ip .eq. 3) THEN !(krokecker delta_p^2)
        xphij   =  imagu * (SQRT2*eta_B + eta_T/SQRT2)
        factj = real(Factorial(ij-1),dp)        
      ELSE
        xphij = 0.; xphijp1 = 0.; xphijm1 = 0.
        factj = 1
      ENDIF 

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
          IF (xNapp2j .NE. (0.,0.)) THEN
            TNapp2j = moments_e(ip+2,ij,ikr,ikz,updatetlevel) * xNapp2j
          ELSE
            TNapp2j = 0.
          ENDIF          
          ! term propto N_e^{p-2,j}
          IF (xNapm2j .NE. (0.,0.)) THEN
            TNapm2j = moments_e(ip-2,ij,ikr,ikz,updatetlevel) * xNapm2j
          ELSE
            TNapm2j = 0.
          ENDIF          
          ! xterm propto N_e^{p,j+1}
          IF (xNapjp1 .NE. (0.,0.)) THEN
            TNapjp1 = moments_e(ip,ij+1,ikr,ikz,updatetlevel) * xNapp2j
          ELSE
            TNapjp1 = 0.
          ENDIF          
          ! term propto N_e^{p,j-1}
          IF (xNapjm1 .NE. (0.,0.)) THEN
            TNapjm1 = moments_e(ip,ij-1,ikr,ikz,updatetlevel) * xNapp2j
          ELSE
            TNapjm1 = 0.
          ENDIF

          ! Collision term completed (Lenard-Bernstein)
          IF (xNapp2j .NE. (0.,0.)) THEN
            TColl = (xColl - nu * b_e2/4.) * moments_e(ip+2,ij,ikr,ikz,updatetlevel)
          ELSE
            TColl = 0.
          ENDIF

          !! Electrical potential term
          Tphi = 0
          IF ( (ip .eq. 1) .or. (ip .eq. 3) ) THEN ! 0 otherwise (krokecker delta_p^0)
            kernelj    = b_e2**(ij-1) * exp(-b_e2)/factj
            kerneljp1  = kernelj * b_e2  /(ij_dp + 1.)
            kerneljm1  = kernelj * ij_dp / b_e2
            Tphi = (xphij * Kernelj   + xphijp1 * Kerneljp1 + xphijm1 * Kerneljm1)* kz * phi(ikr,ikz)
          ENDIF

          ! Sum of all precomputed terms
          moments_rhs_e(ip,ij,ikr,ikz,updatetlevel) = &
              kz * (TNapj + TNapp2j + TNapm2j + TNapjp1 + TNapjm1) &
              + Tphi + TColl

        END DO kzloope
      END DO krloope

    END DO jloope
  END DO ploope


  !!!!!!!!! Ions moments RHS !!!!!!!!!
  Tphi    = 0 ! electrostatic potential term
  taui_qi_etaBi = tau_i/real(q_i)*eta_B * imagu ! one-time computable factor
  xb_i2 = sigma_i**2 * tau_i/2.0 ! species dependant factor of the Kernel, squared

  ploopi : DO ip = ips_i, ipe_i
    ip_dp = real(ip-1,dp) 

    ! x N_i^{p+2,j}
    IF (ip+2 .LE. pmaxi + 1) THEN
      xNapp2j = taui_qi_etaBi * SQRT(ip_dp + 1.) * SQRT(ip_dp + 2.)
    ELSE
      xNapp2j = 0.
    ENDIF

    ! x N_i^{p-2,j}
    IF (ip-2 .GE. 1) THEN
      xNapm2j = taui_qi_etaBi * SQRT(ip_dp) * SQRT(ip_dp - 1.)
    ELSE
      xNapm2j = 0.
    ENDIF

    jloopi : DO ij = ijs_i, ije_i
      ij_dp = real(ij-1,dp) 

      ! x N_i^{p,j+1}
      IF (ij+1 .LE. jmaxi+1) THEN
        xNapjp1 = -taui_qi_etaBi * (ij_dp + 1.)
      ELSE
        xNapjp1 = 0.
      ENDIF

      ! x N_i^{p,j-1}
      IF (ij-1 .GE. 1) THEN
        xNapjm1 = -taui_qi_etaBi * ij_dp
      ELSE
        xNapjm1 = 0.
      ENDIF

      ! x N_i^{pj}
      xNapj   = taui_qi_etaBi * 2.*(ip_dp + ij_dp + 1.)

      ! first part of collision operator (Lenard-Bernstein)
      xColl = -nu*(ip_dp + 2.0*ij_dp)

      ! x phi
      IF (ip .EQ. 1) THEN !(krokecker delta_p^0)
        xphij   =  imagu * (eta_n + 2.*ij_dp*eta_T + 2.*eta_B*(ij_dp+1.) )
        xphijp1 = -imagu * (eta_B + eta_T)*(ij_dp+1.)
        xphijm1 = -imagu * (eta_B + eta_T)* ij_dp
        factj = REAL(Factorial(ij-1),dp)
      ELSE IF (ip .EQ. 3) THEN !(krokecker delta_p^2)
        xphij   = imagu * (SQRT2*eta_B + eta_T/SQRT2)
        factj   = REAL(Factorial(ij-1),dp)        
      ELSE
        xphij = 0.; xphijp1 = 0.; xphijm1 = 0.
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
          
          IF (xNapp2j .NE. (0.,0.)) THEN ! term propto N_i^{p+2,j}
            TNapp2j = moments_i(ip+2,ij,ikr,ikz,updatetlevel) * xNapp2j
          ELSE
            TNapp2j = 0.
          ENDIF          
          IF (xNapm2j .NE. (0.,0.)) THEN ! term propto N_i^{p-2,j}
            TNapm2j = moments_i(ip-2,ij,ikr,ikz,updatetlevel) * xNapm2j
          ELSE
            TNapm2j = 0.
          ENDIF          
          IF (xNapjp1 .NE. (0.,0.)) THEN ! xterm propto N_a^{p,j+1}
            TNapjp1 = moments_i(ip,ij+1,ikr,ikz,updatetlevel) * xNapp2j
          ELSE
            TNapjp1 = 0.
          ENDIF          
          IF (xNapjm1 .NE. (0.,0.)) THEN ! term propto N_a^{p,j-1}
            TNapjm1 = moments_i(ip,ij-1,ikr,ikz,updatetlevel) * xNapp2j
          ELSE
            TNapjm1 = 0.
          ENDIF

          ! Collision term completed (Lenard-Bernstein)
          IF (xNapp2j .NE. (0.,0.)) THEN
            TColl = (xColl - nu * b_i2/4.) * moments_i(ip+2,ij,ikr,ikz,updatetlevel)
          ELSE
            TColl = 0.
          ENDIF

          !! Electrical potential term
          Tphi = 0
          IF ( (ip .eq. 1) .or. (ip .eq. 3) ) THEN ! 0 otherwise (krokecker delta_p^0, delta_p^2)
            kernelj    = b_i2**(ij-1) * exp(-b_i2)/factj
            kerneljp1  = kernelj * b_i2  /(ij_dp + 1._dp)
            kerneljm1  = kernelj * ij_dp / b_i2
            Tphi = (xphij * Kernelj   + xphijp1 * Kerneljp1 + xphijm1 * Kerneljm1)* kz * phi(ikr,ikz)
          ENDIF

          ! Sum of all precomputed terms
          moments_rhs_i(ip,ij,ikr,ikz,updatetlevel) = &
              kz * (TNapj + TNapp2j + TNapm2j + TNapjp1 + TNapjm1) &
              + Tphi + TColl

        END DO kzloopi
      END DO krloopi

    END DO jloopi
  END DO ploopi

END SUBROUTINE moments_eq_rhs
