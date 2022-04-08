MODULE moments_eq_rhs
  IMPLICIT NONE
  PUBLIC :: compute_moments_eq_rhs
CONTAINS

SUBROUTINE compute_moments_eq_rhs
  USE model, only: KIN_E
  IMPLICIT NONE
  IF(KIN_E) CALL moments_eq_rhs_e
            CALL moments_eq_rhs_i
END SUBROUTINE compute_moments_eq_rhs
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
  USE calculus, ONLY : interp_z, grad_z, grad_z2
  IMPLICIT NONE

  INTEGER     :: p_int, j_int ! loops indices and polynom. degrees
  REAL(dp)    :: kx, ky, kperp2, dzlnB_o_J
  COMPLEX(dp) :: Tnepj, Tnepp2j, Tnepm2j, Tnepjp1, Tnepjm1 ! Terms from b x gradB and drives
  COMPLEX(dp) :: Tnepp1j, Tnepm1j, Tnepp1jm1, Tnepm1jm1 ! Terms from mirror force with non adiab moments
  COMPLEX(dp) :: Tperp, Tpar, Tmir, Tphi
  COMPLEX(dp) :: Unepm1j, Unepm1jp1, Unepm1jm1 ! Terms from mirror force with adiab moments
  COMPLEX(dp) :: i_ky
  ! To store derivatives and odd-even z grid interpolations
  COMPLEX(dp), DIMENSION(izs:ize) :: ddznepp1j, ddznepm1j, &
              nepp1j, nepp1jm1, nepm1j, nepm1jm1, nepm1jp1, ddz2Nepj

   ! Measuring execution time
  CALL cpu_time(t0_rhs)

  ! Kinetic loops
  ploope : DO ip = ips_e, ipe_e ! loop over Hermite degree indices
    p_int = parray_e(ip)    ! Hermite polynom. degree
    eo    = MODULO(p_int,2) ! Indicates if we are on even or odd z grid
    jloope : DO ij = ijs_e, ije_e ! loop over Laguerre degree indices
    j_int = jarray_e(ij)
      ! Spatial loops
      kxloope : DO ikx = ikxs,ikxe
      kx     = kxarray(ikx)   ! radial wavevector
      kyloope : DO iky = ikys,ikye
      ky     = kyarray(iky)   ! toroidal wavevector
        i_ky   = imagu * ky     ! toroidal derivative
        IF (Nky .EQ. 1) i_ky = imagu * kxarray(ikx) ! If 1D simulation we put kx as ky
        ! Compute z derivatives and odd-even z interpolations
        CALL   grad_z(eo,nadiab_moments_e(ip+1,ij  ,ikx,iky,izgs:izge), ddznepp1j(izs:ize))
        CALL   grad_z(eo,nadiab_moments_e(ip-1,ij  ,ikx,iky,izgs:izge), ddznepm1j(izs:ize))
        CALL interp_z(eo,nadiab_moments_e(ip+1,ij  ,ikx,iky,izgs:izge),  nepp1j  (izs:ize))
        CALL interp_z(eo,nadiab_moments_e(ip+1,ij-1,ikx,iky,izgs:izge),  nepp1jm1(izs:ize))
        CALL interp_z(eo,nadiab_moments_e(ip-1,ij  ,ikx,iky,izgs:izge),  nepm1j  (izs:ize))
        CALL interp_z(eo,nadiab_moments_e(ip-1,ij+1,ikx,iky,izgs:izge),  nepm1jp1(izs:ize))
        CALL interp_z(eo,nadiab_moments_e(ip-1,ij-1,ikx,iky,izgs:izge),  nepm1jm1(izs:ize))
        ! Parallel hyperdiffusion
        CALL  grad_z2(moments_e(ip,ij,ikx,iky,izgs:izge,updatetlevel),   ddz2Nepj(izs:ize))

        zloope : DO  iz = izs,ize
          ! kperp
          kperp2= kparray(ikx,iky,iz,eo)**2

          !! Compute moments mixing terms
          Tperp = 0._dp; Tpar = 0._dp; Tmir = 0._dp
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
          ! ddz derivative for Landau damping term
          Tpar = xnepp1j(ip) * ddznepp1j(iz) + xnepm1j(ip) * ddznepm1j(iz)
          ! Mirror terms (trapping)
          Tnepp1j   = ynepp1j  (ip,ij) * nepp1j  (iz)
          Tnepp1jm1 = ynepp1jm1(ip,ij) * nepp1jm1(iz)
          Tnepm1j   = ynepm1j  (ip,ij) * nepm1j  (iz)
          Tnepm1jm1 = ynepm1jm1(ip,ij) * nepm1jm1(iz)
          ! Trapping terms
          Unepm1j   = znepm1j  (ip,ij) * nepm1j  (iz)
          Unepm1jp1 = znepm1jp1(ip,ij) * nepm1jp1(iz)
          Unepm1jm1 = znepm1jm1(ip,ij) * nepm1jm1(iz)

          Tmir = Tnepp1j + Tnepp1jm1 + Tnepm1j + Tnepm1jm1 + Unepm1j + Unepm1jp1 + Unepm1jm1
          !! Electrical potential term
          IF ( p_int .LE. 2 ) THEN ! kronecker p0 p1 p2
            Tphi = (xphij_i  (ip,ij)*kernel_e(ij  ,ikx,iky,iz,eo) &
                  + xphijp1_i(ip,ij)*kernel_e(ij+1,ikx,iky,iz,eo) &
                  + xphijm1_i(ip,ij)*kernel_e(ij-1,ikx,iky,iz,eo))*phi(ikx,iky,iz)
          ELSE
            Tphi = 0._dp
          ENDIF

          !! Sum of all RHS terms
          moments_rhs_e(ip,ij,ikx,iky,iz,updatetlevel) = &
              ! Perpendicular magnetic gradient/curvature effects
              - imagu*Ckxky(ikx,iky,iz,eo)*hatR(iz,eo)* (Tnepj + Tnepp2j + Tnepm2j + Tnepjp1 + Tnepjm1)&
              ! Parallel coupling (Landau Damping)
              - Tpar*gradz_coeff(iz,eo) &
              ! Mirror term (parallel magnetic gradient)
              - gradzB(iz,eo)* Tmir  *gradz_coeff(iz,eo) &
              ! Drives (density + temperature gradients)
              - i_ky * Tphi &
              ! Numerical perpendicular hyperdiffusion (totally artificial, for stability purpose)
              - (mu_x*kx**4 + mu_y*ky**4)*moments_e(ip,ij,ikx,iky,iz,updatetlevel) &
              ! Numerical parallel hyperdiffusion "+ (mu_z*kz**4)"
              + mu_z * diff_dz_coeff * ddz2Nepj(iz) &
              ! Collision term
              + TColl_e(ip,ij,ikx,iky,iz) &
              ! Nonlinear term
              - Sepj(ip,ij,ikx,iky,iz)

         END DO zloope
        END DO kyloope
      END DO kxloope
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
  USE calculus, ONLY : interp_z, grad_z, grad_z2
  IMPLICIT NONE

  INTEGER     :: p_int, j_int ! loops indices and polynom. degrees
  REAL(dp)    :: kx, ky, kperp2
  COMPLEX(dp) :: Tnipj, Tnipp2j, Tnipm2j, Tnipjp1, Tnipjm1
  COMPLEX(dp) :: Tnipp1j, Tnipm1j, Tnipp1jm1, Tnipm1jm1 ! Terms from mirror force with non adiab moments
  COMPLEX(dp) :: Unipm1j, Unipm1jp1, Unipm1jm1 ! Terms from mirror force with adiab moments
  COMPLEX(dp) :: Tperp, Tpar, Tmir, Tphi
  COMPLEX(dp) :: i_ky
  ! To store derivatives and odd-even z grid interpolations
  COMPLEX(dp), DIMENSION(izs:ize) :: ddznipp1j, ddznipm1j, &
              nipp1j, nipp1jm1, nipm1j, nipm1jm1, nipm1jp1, ddz2Nipj
  ! Measuring execution time
  CALL cpu_time(t0_rhs)

  ! Kinetic loops
  ploopi : DO ip = ips_i, ipe_i  ! Hermite loop
    p_int = parray_i(ip)    ! Hermite degree
    eo    = MODULO(p_int,2) ! Indicates if we are on odd or even z grid
    jloopi : DO ij = ijs_i, ije_i  ! This loop is from 1 to jmaxi+1
      j_int = jarray_i(ij)
      ! Spatial loops
      kxloopi : DO ikx = ikxs,ikxe
      kx     = kxarray(ikx)   ! radial wavevector
      kyloopi : DO iky = ikys,ikye
      ky     = kyarray(iky)   ! toroidal wavevector
        i_ky   = imagu * ky     ! toroidal derivative
        IF (Nky .EQ. 1) i_ky = imagu * kxarray(ikx) ! If 1D simulation we put kx as ky
        ! Compute z derivatives and odd-even z interpolations
        CALL   grad_z(eo,nadiab_moments_i(ip+1,ij  ,ikx,iky,izgs:izge),ddznipp1j(izs:ize))
        CALL   grad_z(eo,nadiab_moments_i(ip-1,ij  ,ikx,iky,izgs:izge),ddznipm1j(izs:ize))
        CALL interp_z(eo,nadiab_moments_i(ip+1,ij  ,ikx,iky,izgs:izge), nipp1j  (izs:ize))
        CALL interp_z(eo,nadiab_moments_i(ip+1,ij-1,ikx,iky,izgs:izge), nipp1jm1(izs:ize))
        CALL interp_z(eo,nadiab_moments_i(ip-1,ij  ,ikx,iky,izgs:izge), nipm1j  (izs:ize))
        CALL interp_z(eo,nadiab_moments_i(ip-1,ij+1,ikx,iky,izgs:izge), nipm1jp1(izs:ize))
        CALL interp_z(eo,nadiab_moments_i(ip-1,ij-1,ikx,iky,izgs:izge), nipm1jm1(izs:ize))
        ! for hyperdiffusion
        CALL  grad_z2(moments_i(ip,ij,ikx,iky,izgs:izge,updatetlevel),  ddz2Nipj(:))

        zloopi : DO  iz = izs,ize
          kperp2= kparray(ikx,iky,iz,eo)**2

          !! Compute moments mixing terms
          Tperp = 0._dp; Tpar = 0._dp; Tmir = 0._dp
          ! Perpendicular dynamic
          ! term propto n_i^{p,j}
          Tnipj   = xnipj(ip,ij) * nadiab_moments_i(ip    ,ij  ,ikx,iky,iz)
          ! term propto n_i^{p+2,j}
          Tnipp2j = xnipp2j(ip)  * nadiab_moments_i(ip+pp2,ij  ,ikx,iky,iz)
          ! term propto n_i^{p-2,j}
          Tnipm2j = xnipm2j(ip)  * nadiab_moments_i(ip-pp2,ij  ,ikx,iky,iz)
          ! term propto n_i^{p,j+1}
          Tnipjp1 = xnipjp1(ij)  * nadiab_moments_i(ip    ,ij+1,ikx,iky,iz)
          ! term propto n_i^{p,j-1}
          Tnipjm1 = xnipjm1(ij)  * nadiab_moments_i(ip    ,ij-1,ikx,iky,iz)
          ! Tperp
          Tperp = Tnipj + Tnipp2j + Tnipm2j + Tnipjp1 + Tnipjm1
          ! Parallel dynamic
          ! ddz derivative for Landau damping term
          Tpar = xnipp1j(ip) * ddznipp1j(iz) + xnipm1j(ip) * ddznipm1j(iz)
          ! Mirror terms
          Tnipp1j   = ynipp1j  (ip,ij) * nipp1j  (iz)
          Tnipp1jm1 = ynipp1jm1(ip,ij) * nipp1jm1(iz)
          Tnipm1j   = ynipm1j  (ip,ij) * nipm1j  (iz)
          Tnipm1jm1 = ynipm1jm1(ip,ij) * nipm1jm1(iz)
          ! Trapping terms
          Unipm1j   = znipm1j  (ip,ij) * nipm1j  (iz)
          Unipm1jp1 = znipm1jp1(ip,ij) * nipm1jp1(iz)
          Unipm1jm1 = znipm1jm1(ip,ij) * nipm1jm1(iz)

          Tmir = Tnipp1j + Tnipp1jm1 + Tnipm1j + Tnipm1jm1 + Unipm1j + Unipm1jp1 + Unipm1jm1

          !! Electrical potential term
          IF ( p_int .LE. 2 ) THEN ! kronecker p0 p1 p2
            Tphi = (xphij_i  (ip,ij)*kernel_i(ij  ,ikx,iky,iz,eo) &
                  + xphijp1_i(ip,ij)*kernel_i(ij+1,ikx,iky,iz,eo) &
                  + xphijm1_i(ip,ij)*kernel_i(ij-1,ikx,iky,iz,eo))*phi(ikx,iky,iz)
          ELSE
            Tphi = 0._dp
          ENDIF

          !! Sum of all RHS terms
          moments_rhs_i(ip,ij,ikx,iky,iz,updatetlevel) = &
              ! Perpendicular magnetic gradient/curvature effects
              - imagu*Ckxky(ikx,iky,iz,eo)*hatR(iz,eo) * Tperp &
              ! Parallel coupling (Landau damping)
              - gradz_coeff(iz,eo) * Tpar &
              ! Mirror term (parallel magnetic gradient)
              - gradzB(iz,eo) * gradz_coeff(iz,eo) * Tmir &
              ! Drives (density + temperature gradients)
              - i_ky * Tphi &
              ! Numerical hyperdiffusion (totally artificial, for stability purpose)
              - (mu_x*kx**4 + mu_y*ky**4)*moments_i(ip,ij,ikx,iky,iz,updatetlevel) &
              ! Numerical parallel hyperdiffusion "+ (mu_z*kz**4)"
              + mu_z * diff_dz_coeff * ddz2Nipj(iz) &
              ! Collision term
              + TColl_i(ip,ij,ikx,iky,iz)&
              ! Nonlinear term
              - Sipj(ip,ij,ikx,iky,iz)

          END DO zloopi
        END DO kyloopi
      END DO kxloopi
    END DO jloopi
  END DO ploopi

  ! Execution time end
  CALL cpu_time(t1_rhs)
  tc_rhs = tc_rhs + (t1_rhs-t0_rhs)

END SUBROUTINE moments_eq_rhs_i

END MODULE moments_eq_rhs
