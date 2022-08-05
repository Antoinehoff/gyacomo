MODULE moments_eq_rhs
  IMPLICIT NONE
  PUBLIC :: compute_moments_eq_rhs
CONTAINS

SUBROUTINE compute_moments_eq_rhs
  USE model, only: KIN_E
  IMPLICIT NONE
  IF(KIN_E) CALL moments_eq_rhs_e
            CALL moments_eq_rhs_i
  !! BETA TESTING !!
  IF(.false.) CALL add_Maxwellian_background_terms
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
  COMPLEX(dp) :: Tperp, Tpar, Tmir, Tphi, Tpsi
  COMPLEX(dp) :: Unepm1j, Unepm1jp1, Unepm1jm1 ! Terms from mirror force with adiab moments
  COMPLEX(dp) :: i_kx,i_ky,phikykxz, psikykxz

   ! Measuring execution time
  CALL cpu_time(t0_rhs)

  ! Spatial loops
  zloope : DO  iz = izs,ize
    kxloope : DO ikx = ikxs,ikxe
      kx       = kxarray(ikx)   ! radial wavevector
      i_kx     = imagu * kx     ! toroidal derivative

      kyloope : DO iky = ikys,ikye
        ky     = kyarray(iky)   ! toroidal wavevector
        i_ky   = imagu * ky     ! toroidal derivative
        phikykxz = phi(iky,ikx,iz)! tmp phi value
        psikykxz = psi(iky,ikx,iz)! tmp psi value

        ! Kinetic loops
        jloope : DO ij = ijs_e, ije_e  ! This loop is from 1 to jmaxi+1
          j_int = jarray_e(ij)

          ploope : DO ip = ips_e, ipe_e  ! Hermite loop
            p_int = parray_e(ip)    ! Hermite degree
            eo    = MODULO(p_int,2) ! Indicates if we are on odd or even z grid
            kperp2= kparray(iky,ikx,iz,eo)**2

          IF((CLOS .NE. 1) .OR. (p_int+2*j_int .LE. dmaxe)) THEN
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
              Tphi = (xphij_e  (ip,ij)*kernel_e(ij  ,iky,ikx,iz,eo) &
                    + xphijp1_e(ip,ij)*kernel_e(ij+1,iky,ikx,iz,eo) &
                    + xphijm1_e(ip,ij)*kernel_e(ij-1,iky,ikx,iz,eo))*phikykxz
            ELSE
              Tphi = 0._dp
            ENDIF

            !! Vector potential term
            IF ( (p_int .LE. 3) .AND. (p_int .GT. 1) ) THEN ! Kronecker p1 or p3
              Tpsi = (xpsij_e  (ip,ij)*kernel_e(ij  ,iky,ikx,iz,eo) &
                    + xpsijp1_e(ip,ij)*kernel_e(ij+1,iky,ikx,iz,eo) &
                    + xpsijm1_e(ip,ij)*kernel_e(ij-1,iky,ikx,iz,eo))*psikykxz
            ELSE
              Tpsi = 0._dp
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
                - i_ky * (Tphi - Tpsi) &
                ! Numerical perpendicular hyperdiffusion (totally artificial, for stability purpose)
                - mu_x*diff_kx_coeff*kx**N_HD*moments_e(ip,ij,iky,ikx,iz,updatetlevel) &
                - mu_y*diff_ky_coeff*ky**N_HD*moments_e(ip,ij,iky,ikx,iz,updatetlevel) &
                ! Numerical parallel hyperdiffusion "+ (mu_z*kz**4)"  see Pueschel 2010 (eq 25)
                - mu_z * diff_dz_coeff * ddz4_Nepj(ip,ij,iky,ikx,iz) &
                ! Collision term
                + TColl_e(ip,ij,iky,ikx,iz) &
                ! Nonlinear term
                - hatB_NL(iz,eo) * Sepj(ip,ij,iky,ikx,iz)

            IF(ip-4 .GT. 0) &
              ! Numerical parallel velocity hyperdiffusion "+ dvpar4 g_a" see Pueschel 2010 (eq 33)
              moments_rhs_e(ip,ij,iky,ikx,iz,updatetlevel) = &
                moments_rhs_e(ip,ij,iky,ikx,iz,updatetlevel) &
                  + mu_p * moments_e(ip-4,ij,iky,ikx,iz,updatetlevel)

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
  COMPLEX(dp) :: Tperp, Tpar, Tmir, Tphi, Tpsi
  COMPLEX(dp) :: i_kx, i_ky, phikykxz, psikykxz

  ! Measuring execution time
  CALL cpu_time(t0_rhs)


  ! Spatial loops
  zloopi : DO  iz = izs,ize
    kxloopi : DO ikx = ikxs,ikxe
      kx       = kxarray(ikx)   ! radial wavevector
      i_kx     = imagu * kx     ! toroidal derivative
        kyloopi : DO iky = ikys,ikye
          ky     = kyarray(iky)   ! toroidal wavevector
          i_ky   = imagu * ky     ! toroidal derivative
          phikykxz = phi(iky,ikx,iz)! tmp phi value
          psikykxz = psi(iky,ikx,iz)! tmp phi value

        ! Kinetic loops
        jloopi : DO ij = ijs_i, ije_i  ! This loop is from 1 to jmaxi+1
          j_int = jarray_i(ij)

          ploopi : DO ip = ips_i, ipe_i  ! Hermite loop
            p_int = parray_i(ip)    ! Hermite degree
            eo    = MODULO(p_int,2) ! Indicates if we are on odd or even z grid
            kperp2= kparray(iky,ikx,iz,eo)**2
            IF((CLOS .NE. 1) .OR. (p_int+2*j_int .LE. dmaxi)) THEN
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

              !! Vector potential term
              IF ( (p_int .LE. 3) .AND. (p_int .GT. 1) ) THEN ! Kronecker p1 or p3
                Tpsi = (xpsij_i  (ip,ij)*kernel_i(ij  ,iky,ikx,iz,eo) &
                      + xpsijp1_i(ip,ij)*kernel_i(ij+1,iky,ikx,iz,eo) &
                      + xpsijm1_i(ip,ij)*kernel_i(ij-1,iky,ikx,iz,eo))*psikykxz
              ELSE
                Tpsi = 0._dp
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
                  - i_ky * (Tphi - Tpsi) &
                  ! Numerical hyperdiffusion (totally artificial, for stability purpose)
                  - mu_x*diff_kx_coeff*kx**N_HD*moments_i(ip,ij,iky,ikx,iz,updatetlevel) &
                  - mu_y*diff_ky_coeff*ky**N_HD*moments_i(ip,ij,iky,ikx,iz,updatetlevel) &
                  ! Numerical parallel hyperdiffusion "+ (mu_z*kz**4)"
                  + mu_z * diff_dz_coeff * ddz4_Nipj(ip,ij,iky,ikx,iz) &
                  ! Collision term
                  + TColl_i(ip,ij,iky,ikx,iz)&
                  ! Nonlinear term with a (gxx*gxy - gxy**2)^1/2 factor
                  - hatB_NL(iz,eo) * Sipj(ip,ij,iky,ikx,iz)

                  IF(ip-4 .GT. 0) &
                    ! Numerical parallel velocity hyperdiffusion "+ dvpar4 g_a" see Pueschel 2010 (eq 33)
                    moments_rhs_i(ip,ij,iky,ikx,iz,updatetlevel) = &
                      moments_rhs_i(ip,ij,iky,ikx,iz,updatetlevel) &
                        + mu_p * moments_i(ip-4,ij,iky,ikx,iz,updatetlevel)
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

SUBROUTINE add_Maxwellian_background_terms
  ! This routine is meant to add the terms rising from the magnetic operator,
  ! i.e. (B x GradB) Grad, applied on the background Maxwellian distribution
  ! (x_a + spar^2)(b x GradB) GradFaM
  ! It gives birth to kx=ky=0 sources terms (averages) that hit moments 00, 20,
  ! 40, 01,02, 21 with background gradient dependences.
  USE prec_const
  USE time_integration, ONLY : updatetlevel
  USE model,      ONLY: taue_qe, taui_qi, k_N, eta_N, k_T, eta_T, KIN_E
  USE array,      ONLY: moments_rhs_e, moments_rhs_i
  USE grid,       ONLY: contains_kx0, contains_ky0, ikx_0, iky_0,&
                        ips_e,ipe_e,ijs_e,ije_e,ips_i,ipe_i,ijs_i,ije_i,&
                        jarray_e,parray_e,jarray_i,parray_i, zarray, izs,ize,&
                        ip,ij
  IMPLICIT NONE
  real(dp), DIMENSION(izs:ize) :: sinz

  sinz(izs:ize) = SIN(zarray(izs:ize,0))

  IF(contains_kx0 .AND. contains_ky0) THEN
    IF(KIN_E) THEN
      DO ip = ips_e,ipe_e
        DO ij = ijs_e,ije_e
          SELECT CASE(ij-1)
          CASE(0) ! j = 0
            SELECT CASE (ip-1)
            CASE(0) ! Na00 term
                moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                  +taue_qe * sinz(izs:ize) * (1.5_dp*eta_N*k_N - 1.125_dp*eta_N*k_T)
            CASE(2) ! Na20 term
                moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                  +taue_qe * sinz(izs:ize) * (SQRT2*0.5_dp*eta_N*k_N - 2.75_dp*eta_T*k_T)
            CASE(4) ! Na40 term
                moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                  +taue_qe * sinz(izs:ize) * SQRT6*0.75_dp*eta_T*k_T
            END SELECT
          CASE(1) ! j = 1
            SELECT CASE (ip-1)
            CASE(0) ! Na01 term
                moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                  -taue_qe * sinz(izs:ize) * (eta_N*k_N + 3.5_dp*eta_T*k_T)
            CASE(2) ! Na21 term
                moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                  -taue_qe * sinz(izs:ize) * SQRT2*eta_T*k_T
            END SELECT
          CASE(2) ! j = 2
            SELECT CASE (ip-1)
            CASE(0) ! Na02 term
                moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                  +taue_qe * sinz(izs:ize) * 2._dp*eta_T*k_T
            END SELECT
          END SELECT
        ENDDO
      ENDDO
    ENDIF

    DO ip = ips_i,ipe_i
      DO ij = ijs_i,ije_i
        SELECT CASE(ij-1)
        CASE(0) ! j = 0
          SELECT CASE (ip-1)
          CASE(0) ! Na00 term
              moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                +taui_qi * sinz(izs:ize) * (1.5_dp*k_N - 1.125_dp*k_T)
          CASE(2) ! Na20 term
              moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                +taui_qi * sinz(izs:ize) * (SQRT2*0.5_dp*k_N - 2.75_dp*k_T)
          CASE(4) ! Na40 term
              moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                +taui_qi * sinz(izs:ize) * SQRT6*0.75_dp*k_T
          END SELECT
        CASE(1) ! j = 1
          SELECT CASE (ip-1)
          CASE(0) ! Na01 term
              moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                -taui_qi * sinz(izs:ize) * (k_N + 3.5_dp*k_T)
          CASE(2) ! Na21 term
              moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                -taui_qi * sinz(izs:ize) * SQRT2*k_T
          END SELECT
        CASE(2) ! j = 2
          SELECT CASE (ip-1)
          CASE(0) ! Na02 term
              moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                +taui_qi * sinz(izs:ize) * 2._dp*k_T
          END SELECT
        END SELECT
      ENDDO
    ENDDO

  ENDIF

END SUBROUTINE

END MODULE moments_eq_rhs
