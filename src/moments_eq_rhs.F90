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
  USE collision
  use geometry
  IMPLICIT NONE

  INTEGER     :: p_int, j_int ! loops indices and polynom. degrees
  REAL(dp)    :: kx, ky, kperp2, dzlnB_o_J
  COMPLEX(dp) :: Tnepj, Tnepp2j, Tnepm2j, Tnepjp1, Tnepjm1, Tpare, Tphi ! Terms from b x gradB and drives
  COMPLEX(dp) :: Tmir, Tnepp1j, Tnepm1j, Tnepp1jm1, Tnepm1jm1 ! Terms from mirror force with non adiab moments
  COMPLEX(dp) :: UNepm1j, UNepm1jp1, UNepm1jm1 ! Terms from mirror force with adiab moments
  COMPLEX(dp) :: TColl ! terms of the rhs
  COMPLEX(dp) :: i_ky
  REAL(dp)    :: delta_p0, delta_p1, delta_p2
  INTEGER     :: izm2, izm1, izp1, izp2 ! indices for centered FDF ddz

   ! Measuring execution time
  CALL cpu_time(t0_rhs)

  ploope : DO ip = ips_e, ipe_e ! loop over Hermite degree indices
    p_int = parray_e(ip) ! Hermite polynom. degree
    delta_p0 = 0._dp; delta_p1 = 0._dp; delta_p2 = 0._dp
    IF(p_int .EQ. 0) delta_p0 = 1._dp
    IF(p_int .EQ. 1) delta_p1 = 1._dp
    IF(p_int .EQ. 2) delta_p2 = 1._dp

    jloope : DO ij = ijs_e, ije_e ! loop over Laguerre degree indices
    j_int = jarray_e(ij)
      ! Loop on kspace
      zloope : DO  iz = izs,ize
        ! Obtain the index with an array that accounts for boundary conditions
        !   e.g. : 4 stencil with periodic BC, izarray(Nz+2) = 2, izarray(-1) = Nz-1
        izp1 = izarray(iz+1); izp2 = izarray(iz+2);
        izm1 = izarray(iz-1); izm2 = izarray(iz-2);
        !
        kxloope : DO ikx = ikxs,ikxe
        kx     = kxarray(ikx)   ! radial wavevector
          kyloope : DO iky = ikys,ikye
          ky     = kyarray(iky)   ! toroidal wavevector
          i_ky   = imagu * ky     ! toroidal derivative
          IF (Nky .EQ. 1) i_ky = imagu * kxarray(ikx) ! If 1D simulation we put kx as ky
          ! kperp2= gxx(iz)*kx**2 + 2._dp*gxy(iz)*kx*ky + gyy(iz)*ky**2
          kperp2= kparray(ikx,iky,iz)**2

          !! Compute moments mixing terms
          ! Perpendicular dynamic
          ! term propto n_e^{p,j}
          Tnepj   = xnepj(ip,ij)* nadiab_moments_e(ip,ij,ikx,iky,iz)
          ! term propto n_e^{p+2,j}
          Tnepp2j = xnepp2j(ip) * nadiab_moments_e(ip+pp2,ij,ikx,iky,iz)
          ! term propto n_e^{p-2,j}
          Tnepm2j = xnepm2j(ip) * nadiab_moments_e(ip-pp2,ij,ikx,iky,iz)
          ! term propto n_e^{p,j+1}
          Tnepjp1 = xnepjp1(ij) * nadiab_moments_e(ip,ij+1,ikx,iky,iz)
          ! term propto n_e^{p,j-1}
          Tnepjm1 = xnepjm1(ij) * nadiab_moments_e(ip,ij-1,ikx,iky,iz)
          ! Parallel dynamic
          ! ddz derivative for Landau damping term
          Tpare = xnepp1j(ip) * &
            ( onetwelfth*nadiab_moments_e(ip+1,ij,ikx,iky,izm2)&
             -  twothird*nadiab_moments_e(ip+1,ij,ikx,iky,izm1)&
             +  twothird*nadiab_moments_e(ip+1,ij,ikx,iky,izp1)&
             -onetwelfth*nadiab_moments_e(ip+1,ij,ikx,iky,izp2))&
              +xnepm1j(ip) * &
            ( onetwelfth*nadiab_moments_e(ip-1,ij,ikx,iky,izm2)&
             -  twothird*nadiab_moments_e(ip-1,ij,ikx,iky,izm1)&
             +  twothird*nadiab_moments_e(ip-1,ij,ikx,iky,izp1)&
             -onetwelfth*nadiab_moments_e(ip-1,ij,ikx,iky,izp2))
          ! Mirror terms
          Tnepp1j   =   ynepp1j(ip,ij) * nadiab_moments_e(ip+1,ij  ,ikx,iky,iz)
          Tnepp1jm1 = ynepp1jm1(ip,ij) * nadiab_moments_e(ip+1,ij-1,ikx,iky,iz)
          Tnepm1j   =   ynepm1j(ip,ij) * nadiab_moments_e(ip-1,ij  ,ikx,iky,iz)
          Tnepm1jm1 = ynepm1jm1(ip,ij) * nadiab_moments_e(ip-1,ij-1,ikx,iky,iz)
          ! Trapping terms
          UNepm1j   =   zNepm1j(ip,ij) * nadiab_moments_e(ip-1,ij  ,ikx,iky,iz)
          UNepm1jp1 = zNepm1jp1(ip,ij) * nadiab_moments_e(ip-1,ij+1,ikx,iky,iz)
          UNepm1jm1 = zNepm1jm1(ip,ij) * nadiab_moments_e(ip-1,ij-1,ikx,iky,iz)

          Tmir = Tnepp1j + Tnepp1jm1 + Tnepm1j + Tnepm1jm1 + UNepm1j + UNepm1jp1 + UNepm1jm1

          !! Electrical potential term
          IF ( p_int .LE. 2 ) THEN ! kronecker p0 p1 p2
          Tphi = phi(ikx,iky,iz) * (xphij(ip,ij)*kernel_e(ij,ikx,iky,iz) &
                   + xphijp1(ip,ij)*kernel_e(ij+1,ikx,iky,iz) &
                   + xphijm1(ip,ij)*kernel_e(ij-1,ikx,iky,iz) )
          ELSE
            Tphi = 0._dp
          ENDIF

          !! Collision
          IF (CO .EQ. 0) THEN ! Lenard Bernstein
            CALL LenardBernstein_e(ip,ij,ikx,iky,iz,TColl)
          ELSEIF (CO .EQ. 1) THEN ! GK Dougherty
            CALL DoughertyGK_e(ip,ij,ikx,iky,iz,TColl)
          ELSE ! COSOLver matrix
            TColl = TColl_e(ip,ij,ikx,iky,iz)
          ENDIF

          !! Sum of all linear terms (the sign is inverted to match RHS)
          moments_rhs_e(ip,ij,ikx,iky,iz,updatetlevel) = &
              ! Perpendicular magnetic gradient/curvature effects
              - imagu*Ckxky(ikx,iky,iz)*hatR(iz)* (Tnepj + Tnepp2j + Tnepm2j + Tnepjp1 + Tnepjm1)&
              ! Parallel coupling (Landau Damping)
              - Tpare*inv_deltaz*gradz_coeff(iz) &
              ! Mirror term (parallel magnetic gradient)
              - gradzB(iz)* Tmir  *gradz_coeff(iz) &
              ! Drives (density + temperature gradients)
              - i_ky * Tphi &
              ! Electrostatic background gradients
              - i_ky * K_E * moments_e(ip,ij,ikx,iky,iz,updatetlevel) &
              ! Numerical hyperdiffusion (totally artificial, for stability purpose)
              - mu*kperp2**2 * moments_e(ip,ij,ikx,iky,iz,updatetlevel) &
              ! Collision term
              + TColl

          !! Adding non linearity
          IF ( NON_LIN ) THEN
            moments_rhs_e(ip,ij,ikx,iky,iz,updatetlevel) = &
              moments_rhs_e(ip,ij,ikx,iky,iz,updatetlevel) - Sepj(ip,ij,ikx,iky,iz)
          ENDIF

          END DO kyloope
        END DO kxloope
      END DO zloope
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
  USE collision
  IMPLICIT NONE

  INTEGER     :: p_int, j_int ! loops indices and polynom. degrees
  REAL(dp)    :: kx, ky, kperp2
  COMPLEX(dp) :: Tnipj, Tnipp2j, Tnipm2j, Tnipjp1, Tnipjm1, Tpari, Tphi
  COMPLEX(dp) :: Tmir, Tnipp1j, Tnipm1j, Tnipp1jm1, Tnipm1jm1 ! Terms from mirror force with non adiab moments
  COMPLEX(dp) :: UNipm1j, UNipm1jp1, UNipm1jm1 ! Terms from mirror force with adiab moments
  COMPLEX(dp) :: TColl ! terms of the rhs
  COMPLEX(dp) :: i_ky
  REAL(dp)    :: delta_p0, delta_p1, delta_p2
  INTEGER     :: izm2, izm1, izp1, izp2 ! indices for centered FDF ddz

  ! Measuring execution time
  CALL cpu_time(t0_rhs)

  ploopi : DO ip = ips_i, ipe_i  ! Hermite loop
    p_int= parray_i(ip)   ! Hermite degree

    delta_p0 = 0._dp; delta_p1 = 0._dp; delta_p2 = 0._dp
    IF(p_int .EQ. 0) delta_p0 = 1._dp
    IF(p_int .EQ. 1) delta_p1 = 1._dp
    IF(p_int .EQ. 2) delta_p2 = 1._dp

    jloopi : DO ij = ijs_i, ije_i  ! This loop is from 1 to jmaxi+1
      j_int = jarray_i(ij)
      ! Loop on kspace
      zloopi : DO  iz = izs,ize
        ! Obtain the index with an array that accounts for boundary conditions
        !   e.g. : 4 stencil with periodic BC, izarray(Nz+2) = 2, izarray(-1) = Nz-1
        izp1 = izarray(iz+1); izp2 = izarray(iz+2);
        izm1 = izarray(iz-1); izm2 = izarray(iz-2);
        !
        kxloopi : DO ikx = ikxs,ikxe
        kx     = kxarray(ikx)   ! radial wavevector
          kyloopi : DO iky = ikys,ikye
          ky     = kyarray(iky)   ! toroidal wavevector
          i_ky   = imagu * ky     ! toroidal derivative
          IF (Nky .EQ. 1) i_ky = imagu * kxarray(ikx) ! If 1D simulation we put kx as ky
          ! kperp2= gxx(iz)*kx**2 + 2._dp*gxy(iz)*kx*ky + gyy(iz)*ky**2
          kperp2= kparray(ikx,iky,iz)**2

          !! Compute moments mixing terms
          ! Perpendicular dynamic
          ! term propto n_i^{p,j}
          Tnipj   = xnipj(ip,ij) * nadiab_moments_i(ip,ij,ikx,iky,iz)
          ! term propto n_i^{p+2,j}
          Tnipp2j = xnipp2j(ip)  * nadiab_moments_i(ip+pp2,ij,ikx,iky,iz)
          ! term propto n_i^{p-2,j}
          Tnipm2j = xnipm2j(ip)  * nadiab_moments_i(ip-pp2,ij,ikx,iky,iz)
          ! term propto n_e^{p,j+1}
          Tnipjp1 = xnipjp1(ij)  * nadiab_moments_i(ip,ij+1,ikx,iky,iz)
          ! term propto n_e^{p,j-1}
          Tnipjm1 = xnipjm1(ij)  * nadiab_moments_i(ip,ij-1,ikx,iky,iz)
          ! Parallel dynamic
          ! term propto N_i^{p,j+1}, centered FDF
          Tpari = xnipp1j(ip) * &
            ( onetwelfth*nadiab_moments_i(ip+1,ij,ikx,iky,izm2)&
             -  twothird*nadiab_moments_i(ip+1,ij,ikx,iky,izm1)&
             +  twothird*nadiab_moments_i(ip+1,ij,ikx,iky,izp1)&
             -onetwelfth*nadiab_moments_i(ip+1,ij,ikx,iky,izp2))&
              +xnipm1j(ip) * &
            ( onetwelfth*nadiab_moments_i(ip-1,ij,ikx,iky,izm2)&
             -  twothird*nadiab_moments_i(ip-1,ij,ikx,iky,izm1)&
             +  twothird*nadiab_moments_i(ip-1,ij,ikx,iky,izp1)&
             -onetwelfth*nadiab_moments_i(ip-1,ij,ikx,iky,izp2))
          ! Mirror terms
          Tnipp1j   =   ynipp1j(ip,ij) * nadiab_moments_i(ip+1,ij  ,ikx,iky,iz)
          Tnipp1jm1 = ynipp1jm1(ip,ij) * nadiab_moments_i(ip+1,ij-1,ikx,iky,iz)
          Tnipm1j   =   ynipm1j(ip,ij) * nadiab_moments_i(ip-1,ij  ,ikx,iky,iz)
          Tnipm1jm1 = ynipm1jm1(ip,ij) * nadiab_moments_i(ip-1,ij-1,ikx,iky,iz)
          ! Trapping terms
          Unipm1j   =   znipm1j(ip,ij) * nadiab_moments_i(ip-1,ij  ,ikx,iky,iz)
          Unipm1jp1 = znipm1jp1(ip,ij) * nadiab_moments_i(ip-1,ij+1,ikx,iky,iz)
          Unipm1jm1 = znipm1jm1(ip,ij) * nadiab_moments_i(ip-1,ij-1,ikx,iky,iz)

          Tmir = Tnipp1j + Tnipp1jm1 + Tnipm1j + Tnipm1jm1 + UNipm1j + UNipm1jp1 + UNipm1jm1

          !! Electrical potential term
          IF ( p_int .LE. 2 ) THEN ! kronecker p0 p1 p2
            Tphi = phi(ikx,iky,iz) * (xphij(ip,ij)*kernel_i(ij,ikx,iky,iz) &
                     + xphijp1(ip,ij)*kernel_i(ij+1,ikx,iky,iz) &
                     + xphijm1(ip,ij)*kernel_i(ij-1,ikx,iky,iz) )
          ELSE
            Tphi = 0._dp
          ENDIF

          !! Collision
          IF     (CO .EQ. 0) THEN ! Lenard Bernstein
            CALL LenardBernstein_i(ip,ij,ikx,iky,iz,TColl)
          ELSEIF (CO .EQ. 1) THEN ! GK Dougherty
            CALL DoughertyGK_i(ip,ij,ikx,iky,iz,TColl)
          ELSE! COSOLver matrix (Sugama, Coulomb)
            TColl = TColl_i(ip,ij,ikx,iky,iz)
          ENDIF

          !! Sum of all linear terms (the sign is inverted to match RHS)
          moments_rhs_i(ip,ij,ikx,iky,iz,updatetlevel) = &
              ! Perpendicular magnetic gradient/curvature effects
              - imagu*Ckxky(ikx,iky,iz)*hatR(iz)*(Tnipj + Tnipp2j + Tnipm2j + Tnipjp1 + Tnipjm1)&
              ! Parallel coupling (Landau Damping)
              - Tpari*inv_deltaz*gradz_coeff(iz) &
              ! Mirror term (parallel magnetic gradient)
              - gradzB(iz)*Tmir*gradz_coeff(iz) &
              ! Drives (density + temperature gradients)
              - i_ky * Tphi &
              ! Electrostatic background gradients
              - i_ky * K_E * moments_i(ip,ij,ikx,iky,iz,updatetlevel) &
              ! Numerical hyperdiffusion (totally artificial, for stability purpose)
              - mu*kperp2**2 * moments_i(ip,ij,ikx,iky,iz,updatetlevel) &
              ! Collision term
              + TColl

          !! Adding non linearity
          IF ( NON_LIN ) THEN
           moments_rhs_i(ip,ij,ikx,iky,iz,updatetlevel) = &
             moments_rhs_i(ip,ij,ikx,iky,iz,updatetlevel) - Sipj(ip,ij,ikx,iky,iz)
          ENDIF

          END DO kyloopi
        END DO kxloopi
      END DO zloopi
    END DO jloopi
  END DO ploopi

  ! Execution time end
  CALL cpu_time(t1_rhs)
  tc_rhs = tc_rhs + (t1_rhs-t0_rhs)

END SUBROUTINE moments_eq_rhs_i
