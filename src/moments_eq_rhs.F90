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
  IMPLICIT NONE

  INTEGER     :: p_int, j_int ! loops indices and polynom. degrees
  REAL(dp)    :: kx, ky, kperp2
  COMPLEX(dp) :: Tnepj, Tnepp2j, Tnepm2j, Tnepjp1, Tnepjm1, Tpare, Tphi
  COMPLEX(dp) :: TColl ! terms of the rhs
  COMPLEX(dp) :: i_ky
  REAL(dp)    :: delta_p0, delta_p1, delta_p2
  INTEGER     :: izprev,iznext ! indices of previous and next z slice

   ! Measuring execution time
  CALL cpu_time(t0_rhs)

  ploope : DO ip = ips_e, ipe_e ! loop over Hermite degree indices
    p_int= parray_e(ip) ! Hermite polynom. degree

    delta_p0 = 0._dp; delta_p1 = 0._dp; delta_p2 = 0._dp
    IF(p_int .EQ. 0) delta_p0 = 1._dp
    IF(p_int .EQ. 1) delta_p1 = 1._dp
    IF(p_int .EQ. 2) delta_p2 = 1._dp

    jloope : DO ij = ijs_e, ije_e ! loop over Laguerre degree indices
      ! Loop on kspace
      zloope : DO  iz = izs,ize
        ! Periodic BC for first order derivative
        iznext = iz+1; izprev = iz-1;
        IF(iz .EQ.  1) izprev = Nz
        IF(iz .EQ. Nz) iznext = 1

        kxloope : DO ikx = ikxs,ikxe
        kx     = kxarray(ikx)   ! radial wavevector
          kyloope : DO iky = ikys,ikye
          ky     = kyarray(iky)   ! toroidal wavevector
          i_ky   = imagu * ky     ! toroidal derivative
          IF (Nky .EQ. 1) i_ky = imagu * kxarray(ikx) ! If 1D simulation we put kx as ky
          kperp2 = kx**2 + ky**2  ! perpendicular wavevector

          !! Compute moments mixing terms
          ! Perpendicular dynamic
          ! term propto n_e^{p,j}
          Tnepj   = xnepj(ip,ij) * (moments_e(ip,ij,ikx,iky,iz,updatetlevel) &
                               +kernel_e(ij,ikx,iky)*qe_taue*phi(ikx,iky,iz)*delta_p0)
          ! term propto n_e^{p+2,j}
          Tnepp2j = xnepp2j(ip)  * moments_e(ip+2,ij,ikx,iky,iz,updatetlevel)
          ! term propto n_e^{p-2,j}
          Tnepm2j = xnepm2j(ip)  * (moments_e(ip-2,ij,ikx,iky,iz,updatetlevel) &
                               +kernel_e(ij,ikx,iky)*qe_taue*phi(ikx,iky,iz)*delta_p2)
          ! term propto n_e^{p,j+1}
          Tnepjp1 = xnepjp1(ij)  * (moments_e(ip,ij+1,ikx,iky,iz,updatetlevel) &
                               +kernel_e(ij+1,ikx,iky)*qe_taue*phi(ikx,iky,iz)*delta_p0)
          ! term propto n_e^{p,j-1}
          Tnepjm1 = xnepjm1(ij)  * (moments_e(ip,ij-1,ikx,iky,iz,updatetlevel) &
                               +kernel_e(ij-1,ikx,iky)*qe_taue*phi(ikx,iky,iz)*delta_p0)
          ! Parallel dynamic
          ! term propto n_i^{p+1,j}, centered FDF
          Tpare = xnepp1j(ip) * &
              (moments_e(ip+1,ij,ikx,iky,iznext,updatetlevel)&
              -moments_e(ip+1,ij,ikx,iky,izprev,updatetlevel)) &
                +xnepm1j(ip) * &
              (moments_e(ip-1,ij,ikx,iky,iznext,updatetlevel)+kernel_e(ij,ikx,iky)*qe_taue*phi(ikx,iky,iznext)*delta_p1&
              -moments_e(ip-1,ij,ikx,iky,izprev,updatetlevel)-kernel_e(ij,ikx,iky)*qe_taue*phi(ikx,iky,izprev)*delta_p1)

          !! Electrical potential term
          IF ( p_int .LE. 2 ) THEN ! kronecker p0 p1 p2
          Tphi = phi(ikx,iky,iz) * (xphij(ip,ij)*kernel_e(ij, ikx, iky) &
                   + xphijp1(ip,ij)*kernel_e(ij+1, ikx, iky) &
                   + xphijm1(ip,ij)*kernel_e(ij-1, ikx, iky) )
          ELSE
            Tphi = 0._dp
          ENDIF

          !! Collision
          IF (CO .EQ. 0) THEN ! Lenhard Bernstein
            TColl = -nu_ee*(ip+2*ij-3)*moments_e(ip,ij,ikx,iky,iz,updatetlevel)
          ELSEIF (CO .EQ. 1) THEN ! GK Dougherty
            CALL DoughertyGK_e(ip,ij,ikx,iky,iz,TColl)
          ELSE ! COSOLver matrix
            TColl = TColl_e(ip,ij,ikx,iky,iz)
          ENDIF

          !! Sum of all linear terms (the sign is inverted to match RHS)
          moments_rhs_e(ip,ij,ikx,iky,iz,updatetlevel) = &
              - Ckxky(ikx,iky,iz) * (Tnepj + Tnepp2j + Tnepm2j + Tnepjp1 + Tnepjm1)&
              - Tpare/2._dp/deltaz/q0 &
              + i_ky  * Tphi &
              - mu*kperp2**2 * moments_e(ip,ij,ikx,iky,iz,updatetlevel) &
              + TColl

          !! Adding non linearity
          IF ( NON_LIN ) THEN
            moments_rhs_e(ip,ij,ikx,iky,iz,updatetlevel) = &
              moments_rhs_e(ip,ij,ikx,iky,iz,updatetlevel) + Sepj(ip,ij,ikx,iky,iz)
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
  COMPLEX(dp) :: TColl ! terms of the rhs
  COMPLEX(dp) :: i_ky
  REAL(dp)    :: delta_p0, delta_p1, delta_p2
  INTEGER     :: izprev,iznext ! indices of previous and next z slice

  ! Measuring execution time
  CALL cpu_time(t0_rhs)

  ploopi : DO ip = ips_i, ipe_i  ! Hermite loop
    p_int= parray_i(ip)   ! Hermite degree

    delta_p0 = 0._dp; delta_p1 = 0._dp; delta_p2 = 0._dp
    IF(p_int .EQ. 0) delta_p0 = 1._dp
    IF(p_int .EQ. 1) delta_p1 = 1._dp
    IF(p_int .EQ. 2) delta_p2 = 1._dp

    jloopi : DO ij = ijs_i, ije_i  ! This loop is from 1 to jmaxi+1
      ! Loop on kspace
      zloopi : DO  iz = izs,ize
        ! Periodic BC for first order derivative
        iznext = iz+1; izprev = iz-1;
        IF(iz .EQ.  1) izprev = Nz
        IF(iz .EQ. Nz) iznext = 1

        kxloopi : DO ikx = ikxs,ikxe
        kx     = kxarray(ikx)   ! radial wavevector
          kyloopi : DO iky = ikys,ikye
          ky     = kyarray(iky)   ! toroidal wavevector
          i_ky   = imagu * ky     ! toroidal derivative
          IF (Nky .EQ. 1) i_ky = imagu * kxarray(ikx) ! If 1D simulation we put kx as ky
          kperp2 = kx**2 + ky**2  ! perpendicular wavevector

          !! Compute moments mixing terms
          ! Perpendicular dynamic
          ! term propto n_i^{p,j}
          Tnipj   = xnipj(ip,ij)   * (moments_i(ip,ij,ikx,iky,iz,updatetlevel) &
                             +kernel_i(ij,ikx,iky)*qi_taui*phi(ikx,iky,iz)*delta_p0)
          ! term propto n_i^{p+2,j}
          Tnipp2j = xnipp2j(ip) * moments_i(ip+2,ij,ikx,iky,iz,updatetlevel)
          ! term propto n_i^{p-2,j}
          Tnipm2j = xnipm2j(ip) * (moments_i(ip-2,ij,ikx,iky,iz,updatetlevel) &
                             +kernel_i(ij,ikx,iky)*qi_taui*phi(ikx,iky,iz)*delta_p2)
          ! term propto n_e^{p,j+1}
          Tnipjp1 = xnipjp1(ij)  * (moments_i(ip,ij+1,ikx,iky,iz,updatetlevel) &
                               +kernel_i(ij+1,ikx,iky)*qi_taui*phi(ikx,iky,iz)*delta_p0)
          ! term propto n_e^{p,j-1}
          Tnipjm1 = xnipjm1(ij)  * (moments_i(ip,ij-1,ikx,iky,iz,updatetlevel) &
                               +kernel_i(ij-1,ikx,iky)*qi_taui*phi(ikx,iky,iz)*delta_p0)
          ! Parallel dynamic
          ! term propto N_i^{p,j+1}, centered FDF
          Tpari = xnipp1j(ip) * &
              (moments_i(ip+1,ij,ikx,iky,iznext,updatetlevel)&
              -moments_i(ip+1,ij,ikx,iky,izprev,updatetlevel)) &
                +xnipm1j(ip) * &
              (moments_i(ip-1,ij,ikx,iky,iznext,updatetlevel)+qi_taui*kernel_i(ij,ikx,iky)*phi(ikx,iky,iznext)*delta_p1&
              -moments_i(ip-1,ij,ikx,iky,izprev,updatetlevel)-qi_taui*kernel_i(ij,ikx,iky)*phi(ikx,iky,izprev)*delta_p1)

          !! Electrical potential term
          IF ( p_int .LE. 2 ) THEN ! kronecker p0 p1 p2
            Tphi = phi(ikx,iky,iz) * (xphij(ip,ij)*kernel_i(ij, ikx, iky) &
                     + xphijp1(ip,ij)*kernel_i(ij+1, ikx, iky) &
                     + xphijm1(ip,ij)*kernel_i(ij-1, ikx, iky) )
          ELSE
            Tphi = 0._dp
          ENDIF

          !! Collision
          IF     (CO .EQ. 0) THEN ! Lenhard Bernstein
            TColl = -nu_i*(ip+2._dp*ij-3)*moments_i(ip,ij,ikx,iky,iz,updatetlevel)
          ELSEIF (CO .EQ. 1) THEN ! GK Dougherty
            CALL DoughertyGK_i(ip,ij,ikx,iky,iz,TColl)
          ELSE! COSOLver matrix (Sugama, Coulomb)
            TColl = TColl_i(ip,ij,ikx,iky,iz)
          ENDIF

          !! Sum of linear terms
          moments_rhs_i(ip,ij,ikx,iky,iz,updatetlevel) = &
              - Ckxky(ikx,iky,iz) * (Tnipj + Tnipp2j + Tnipm2j + Tnipjp1 + Tnipjm1)&
              - Tpari/2._dp/deltaz/q0 &
              + i_ky  * Tphi &
              - mu*kperp2**2 * moments_i(ip,ij,ikx,iky,iz,updatetlevel) &
              + TColl

          !! Adding non linearity
          IF ( NON_LIN ) THEN
           moments_rhs_i(ip,ij,ikx,iky,iz,updatetlevel) = &
             moments_rhs_i(ip,ij,ikx,iky,iz,updatetlevel) + Sipj(ip,ij,ikx,iky,iz)
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
