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
  COMPLEX(dp) :: i_ky,phikykxz

   ! Measuring execution time
  CALL cpu_time(t0_rhs)

  ! Spatial loops
  zloope : DO  iz = izs,ize
    kxloope : DO ikx = ikxs,ikxe
      kx       = kxarray(ikx)   ! radial wavevector

      kyloope : DO iky = ikys,ikye
        ky     = kyarray(iky)   ! toroidal wavevector
        i_ky   = imagu * ky     ! toroidal derivative
        phikykxz = phi(iky,ikx,iz)! tmp phi value

        ! Kinetic loops
        jloope : DO ij = ijs_e, ije_e  ! This loop is from 1 to jmaxi+1
          j_int = jarray_e(ij)

          ploope : DO ip = ips_e, ipe_e  ! Hermite loop
            p_int = parray_e(ip)    ! Hermite degree
            eo    = MODULO(p_int,2) ! Indicates if we are on odd or even z grid
            kperp2= kparray(iky,ikx,iz,eo)**2

          IF((CLOS .EQ. 1) .OR. (p_int+2*j_int .LE. dmaxe)) THEN
            !! Compute moments mixing terms
            Tperp = 0._dp; Tpar = 0._dp; Tmir = 0._dp
            ! Perpendicular dynamic
            ! term propto n_e^{p,j}
            Tnepj   = xnepj(ip,ij)* nadiab_moments_e(ip,ij,iky,ikx,iz)
            ! term propto n_e^{p+2,j}
            Tnepp2j = xnepp2j(ip) * nadiab_moments_e(ip+pp2,ij,iky,ikx,iz)
            ! term propto n_e^{p-2,j}
            Tnepm2j = xnepm2j(ip) * nadiab_moments_e(ip-pp2,ij,iky,ikx,iz)
            ! term propto n_e^{p,j+1}
            Tnepjp1 = xnepjp1(ij) * nadiab_moments_e(ip,ij+1,iky,ikx,iz)
            ! term propto n_e^{p,j-1}
            Tnepjm1 = xnepjm1(ij) * nadiab_moments_e(ip,ij-1,iky,ikx,iz)
            ! Parallel dynamic
            ! ddz derivative for Landau damping term
            Tpar      = xnepp1j(ip) * ddz_nepj(ip+1,ij,iky,ikx,iz) &
                      + xnepm1j(ip) * ddz_nepj(ip-1,ij,iky,ikx,iz)
            ! Mirror terms
            Tnepp1j   = ynepp1j  (ip,ij) * interp_nepj(ip+1,ij  ,iky,ikx,iz)
            Tnepp1jm1 = ynepp1jm1(ip,ij) * interp_nepj(ip+1,ij-1,iky,ikx,iz)
            Tnepm1j   = ynepm1j  (ip,ij) * interp_nepj(ip-1,ij  ,iky,ikx,iz)
            Tnepm1jm1 = ynepm1jm1(ip,ij) * interp_nepj(ip-1,ij-1,iky,ikx,iz)
            ! Trapping terms
            Unepm1j   = znepm1j  (ip,ij) * interp_nepj(ip-1,ij  ,iky,ikx,iz)
            Unepm1jp1 = znepm1jp1(ip,ij) * interp_nepj(ip-1,ij+1,iky,ikx,iz)
            Unepm1jm1 = znepm1jm1(ip,ij) * interp_nepj(ip-1,ij-1,iky,ikx,iz)

            Tmir = Tnepp1j + Tnepp1jm1 + Tnepm1j + Tnepm1jm1 + Unepm1j + Unepm1jp1 + Unepm1jm1
            !! Electrical potential term
            IF ( p_int .LE. 2 ) THEN ! kronecker p0 p1 p2
              Tphi = (xphij_i  (ip,ij)*kernel_e(ij  ,iky,ikx,iz,eo) &
                    + xphijp1_i(ip,ij)*kernel_e(ij+1,iky,ikx,iz,eo) &
                    + xphijm1_i(ip,ij)*kernel_e(ij-1,iky,ikx,iz,eo))*phikykxz
            ELSE
              Tphi = 0._dp
            ENDIF

            !! Sum of all RHS terms
            moments_rhs_e(ip,ij,iky,ikx,iz,updatetlevel) = &
                ! Perpendicular magnetic gradient/curvature effects
                - imagu*Ckxky(iky,ikx,iz,eo)*hatR(iz,eo)* (Tnepj + Tnepp2j + Tnepm2j + Tnepjp1 + Tnepjm1)&
                ! Parallel coupling (Landau Damping)
                - Tpar*gradz_coeff(iz,eo) &
                ! Mirror term (parallel magnetic gradient)
                - gradzB(iz,eo)* Tmir  *gradz_coeff(iz,eo) &
                ! Drives (density + temperature gradients)
                - i_ky * Tphi &
                ! Numerical perpendicular hyperdiffusion (totally artificial, for stability purpose)
                - (mu_x*kx**4 + mu_y*ky**4)*moments_e(ip,ij,iky,ikx,iz,updatetlevel) &
                ! Numerical parallel hyperdiffusion "+ (mu_z*kz**4)"
                + mu_z * diff_dz_coeff * ddz4_Nepj(ip,ij,iky,ikx,iz) &
                ! Collision term
                + TColl_e(ip,ij,iky,ikx,iz) &
                ! Nonlinear term
                - Sepj(ip,ij,iky,ikx,iz)
          ELSE
            moments_rhs_e(ip,ij,iky,ikx,iz,updatetlevel) = 0._dp
          ENDIF
          END DO ploope
        END DO jloope
      END DO kyloope
    END DO kxloope
  END DO zloope

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
  USE geometry
  USE calculus, ONLY : interp_z, grad_z, grad_z2
  IMPLICIT NONE

  INTEGER     :: p_int, j_int ! loops indices and polynom. degrees
  REAL(dp)    :: kx, ky, kperp2
  COMPLEX(dp) :: Tnipj, Tnipp2j, Tnipm2j, Tnipjp1, Tnipjm1
  COMPLEX(dp) :: Tnipp1j, Tnipm1j, Tnipp1jm1, Tnipm1jm1 ! Terms from mirror force with non adiab moments
  COMPLEX(dp) :: Unipm1j, Unipm1jp1, Unipm1jm1 ! Terms from mirror force with adiab moments
  COMPLEX(dp) :: Tperp, Tpar, Tmir, Tphi
  COMPLEX(dp) :: i_ky, phikykxz

  ! Measuring execution time
  CALL cpu_time(t0_rhs)


  ! Spatial loops
  zloopi : DO  iz = izs,ize
    kxloopi : DO ikx = ikxs,ikxe
      kx       = kxarray(ikx)   ! radial wavevector

        kyloopi : DO iky = ikys,ikye
          ky     = kyarray(iky)   ! toroidal wavevector
          i_ky   = imagu * ky     ! toroidal derivative
          phikykxz = phi(iky,ikx,iz)! tmp phi value

        ! Kinetic loops
        jloopi : DO ij = ijs_i, ije_i  ! This loop is from 1 to jmaxi+1
          j_int = jarray_i(ij)

          ploopi : DO ip = ips_i, ipe_i  ! Hermite loop
            p_int = parray_i(ip)    ! Hermite degree
            eo    = MODULO(p_int,2) ! Indicates if we are on odd or even z grid
            kperp2= kparray(iky,ikx,iz,eo)**2
            IF((CLOS .EQ. 1) .OR. (p_int+2*j_int .LE. dmaxi)) THEN
              !! Compute moments mixing terms
              Tperp = 0._dp; Tpar = 0._dp; Tmir = 0._dp
              ! Perpendicular dynamic
              ! term propto n_i^{p,j}
              Tnipj   = xnipj(ip,ij) * nadiab_moments_i(ip    ,ij  ,iky,ikx,iz)
              ! term propto n_i^{p+2,j}
              Tnipp2j = xnipp2j(ip)  * nadiab_moments_i(ip+pp2,ij  ,iky,ikx,iz)
              ! term propto n_i^{p-2,j}
              Tnipm2j = xnipm2j(ip)  * nadiab_moments_i(ip-pp2,ij  ,iky,ikx,iz)
              ! term propto n_i^{p,j+1}
              Tnipjp1 = xnipjp1(ij)  * nadiab_moments_i(ip    ,ij+1,iky,ikx,iz)
              ! term propto n_i^{p,j-1}
              Tnipjm1 = xnipjm1(ij)  * nadiab_moments_i(ip    ,ij-1,iky,ikx,iz)
              ! Tperp
              Tperp   = Tnipj + Tnipp2j + Tnipm2j + Tnipjp1 + Tnipjm1
              ! Parallel dynamic
              ! ddz derivative for Landau damping term
              Tpar      = xnipp1j(ip) * ddz_nipj(ip+1,ij,iky,ikx,iz) &
                        + xnipm1j(ip) * ddz_nipj(ip-1,ij,iky,ikx,iz)
              ! Mirror terms
              Tnipp1j   = ynipp1j  (ip,ij) * interp_nipj(ip+1,ij  ,iky,ikx,iz)
              Tnipp1jm1 = ynipp1jm1(ip,ij) * interp_nipj(ip+1,ij-1,iky,ikx,iz)
              Tnipm1j   = ynipm1j  (ip,ij) * interp_nipj(ip-1,ij  ,iky,ikx,iz)
              Tnipm1jm1 = ynipm1jm1(ip,ij) * interp_nipj(ip-1,ij-1,iky,ikx,iz)
              ! Trapping terms
              Unipm1j   = znipm1j  (ip,ij) * interp_nipj(ip-1,ij  ,iky,ikx,iz)
              Unipm1jp1 = znipm1jp1(ip,ij) * interp_nipj(ip-1,ij+1,iky,ikx,iz)
              Unipm1jm1 = znipm1jm1(ip,ij) * interp_nipj(ip-1,ij-1,iky,ikx,iz)

              Tmir = Tnipp1j + Tnipp1jm1 + Tnipm1j + Tnipm1jm1 + Unipm1j + Unipm1jp1 + Unipm1jm1

              !! Electrical potential term
              IF ( p_int .LE. 2 ) THEN ! kronecker p0 p1 p2
                Tphi = (xphij_i  (ip,ij)*kernel_i(ij  ,iky,ikx,iz,eo) &
                      + xphijp1_i(ip,ij)*kernel_i(ij+1,iky,ikx,iz,eo) &
                      + xphijm1_i(ip,ij)*kernel_i(ij-1,iky,ikx,iz,eo))*phikykxz
              ELSE
                Tphi = 0._dp
              ENDIF

              !! Sum of all RHS terms
              moments_rhs_i(ip,ij,iky,ikx,iz,updatetlevel) = &
                  ! Perpendicular magnetic gradient/curvature effects
                  - imagu*Ckxky(iky,ikx,iz,eo)*hatR(iz,eo) * Tperp &
                  ! Parallel coupling (Landau damping)
                  - gradz_coeff(iz,eo) * Tpar &
                  ! Mirror term (parallel magnetic gradient)
                  - gradzB(iz,eo) * gradz_coeff(iz,eo) * Tmir &
                  ! Drives (density + temperature gradients)
                  - i_ky * Tphi &
                  ! Numerical hyperdiffusion (totally artificial, for stability purpose)
                  - (mu_x*kx**4 + mu_y*ky**4)*moments_i(ip,ij,iky,ikx,iz,updatetlevel) &
                  ! Numerical parallel hyperdiffusion "+ (mu_z*kz**4)"
                  + mu_z * diff_dz_coeff * ddz4_Nipj(ip,ij,iky,ikx,iz) &
                  ! Collision term
                  + TColl_i(ip,ij,iky,ikx,iz)&
                  ! Nonlinear term
                  - Sipj(ip,ij,iky,ikx,iz)
          ELSE
            moments_rhs_i(ip,ij,iky,ikx,iz,updatetlevel) = 0._dp
          ENDIF
          END DO ploopi
        END DO jloopi
      END DO kyloopi
    END DO kxloopi
  END DO zloopi

  ! Execution time end
  CALL cpu_time(t1_rhs)
  tc_rhs = tc_rhs + (t1_rhs-t0_rhs)

END SUBROUTINE moments_eq_rhs_i

END MODULE moments_eq_rhs
