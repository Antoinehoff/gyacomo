MODULE moments_eq_rhs
  IMPLICIT NONE
  PUBLIC :: compute_moments_eq_rhs
CONTAINS

SUBROUTINE compute_moments_eq_rhs
  USE model
  USE array
  USE fields
  USE grid
  USE basic
  USE prec_const
  USE collision
  USE time_integration
  USE geometry, ONLY: gradz_coeff, dBdz, Ckxky, Gamma_NL!, Gamma_phipar
  USE calculus, ONLY : interp_z, grad_z, grad_z2
  IMPLICIT NONE

    !compute ion moments_eq_rhs
    CALL moments_eq_rhs(ips_i,ipe_i,ipgs_i,ipge_i,ijs_i,ije_i,ijgs_i,ijge_i,jarray_i,parray_i,&
                     xnipj, xnipp2j, xnipm2j, xnipjp1, xnipjm1, xnipp1j, xnipm1j,&
                     ynipp1j, ynipp1jm1, ynipm1j, ynipm1jm1, &
                     znipm1j, znipm1jp1, znipm1jm1, &
                     xphij_i, xphijp1_i, xphijm1_i, xpsij_i, xpsijp1_i, xpsijm1_i,&
                     kernel_i, nadiab_moments_i, ddz_nipj, interp_nipj, Sipj,&
                     moments_i(ipgs_i:ipge_i,ijgs_i:ijge_i,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel),&
                     TColl_i, ddzND_nipj, &
                     moments_rhs_i(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,ikxs:ikxe,izs:ize,updatetlevel))

    !compute ion moments_eq_rhs
    IF(KIN_E) &
    CALL moments_eq_rhs(ips_e,ipe_e,ipgs_e,ipge_e,ijs_e,ije_e,ijgs_e,ijge_e,jarray_e,parray_e,&
                     xnepj, xnepp2j, xnepm2j, xnepjp1, xnepjm1, xnepp1j, xnepm1j,&
                     ynepp1j, ynepp1jm1, ynepm1j, ynepm1jm1, &
                     znepm1j, znepm1jp1, znepm1jm1, &
                     xphij_e, xphijp1_e, xphijm1_e, xpsij_e, xpsijp1_e, xpsijm1_e,&
                     kernel_e, nadiab_moments_e, ddz_nepj, interp_nepj, Sepj,&
                     moments_e(ipgs_e:ipge_e,ijgs_e:ijge_e,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel),&
                     TColl_e, ddzND_nepj,&
                     moments_rhs_e(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,ikxs:ikxe,izs:ize,updatetlevel))

  CONTAINS
  !_____________________________________________________________________________!
  !_____________________________________________________________________________!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! moments_ RHS computation !!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This routine assemble the RHS of the moment hierarchy equations. It uses
  ! linear coefficients that are stored in arrays (xn*, yn* and zn*) computed in
  ! numerics_mod.F90. Otherwise it simply adds the collision term TColl that is
  ! computed in collision_mod.F90 and the nonlinear term Sapj computed in
  ! nonlinear_mod.F90.
  ! All arguments of the subroutines are inputs only except the last one,
  ! moments_rhs_ that will contain the sum of every terms in the RHS.
  !_____________________________________________________________________________!
  SUBROUTINE moments_eq_rhs(ips,ipe,ipgs,ipge,ijs,ije,ijgs,ijge,jarray,parray,&
                   xnapj, xnapp2j, xnapm2j, xnapjp1, xnapjm1, xnapp1j, xnapm1j,&
                   ynapp1j, ynapp1jm1, ynapm1j, ynapm1jm1, &
                   znapm1j, znapm1jp1, znapm1jm1, &
                   xphij, xphijp1, xphijm1, xpsij, xpsijp1, xpsijm1,&
                   kernel, nadiab_moments, ddz_napj, interp_napj, Sapj,&
                   moments_, TColl_, ddzND_napj, moments_rhs_)

    IMPLICIT NONE
    !! INPUTS
    INTEGER, INTENT(IN) :: ips, ipe, ipgs, ipge, ijs, ije, ijgs, ijge
    INTEGER,  DIMENSION(ips:ipe), INTENT(IN) :: parray
    INTEGER,  DIMENSION(ijs:ije), INTENT(IN) :: jarray
    REAL(dp), DIMENSION(ips:ipe,ijs:ije), INTENT(IN) :: xnapj,&
                                        ynapp1j, ynapm1j,   ynapp1jm1, ynapm1jm1,&
                                        znapm1j, znapm1jp1, znapm1jm1,&
                                        xphij, xphijp1, xphijm1,&
                                        xpsij, xpsijp1, xpsijm1
    REAL(dp), DIMENSION(ips:ipe), INTENT(IN) :: &
                                        xnapp1j, xnapm1j,   xnapp2j,   xnapm2j

    REAL(dp), DIMENSION(ijs:ije), INTENT(IN) :: xnapjp1, xnapjm1

    REAL(dp), DIMENSION(ijgs:ijge,ikys:ikye,ikxs:ikxe,izgs:izge,0:1),INTENT(IN) :: kernel

    COMPLEX(dp), DIMENSION(ipgs:ipge,ijgs:ijge,ikys:ikye,ikxs:ikxe,izgs:izge),INTENT(IN) ::&
        moments_,nadiab_moments, ddz_napj, interp_napj, ddzND_napj
    COMPLEX(dp), DIMENSION(ips:ipe,ijs:ije,ikys:ikye,ikxs:ikxe,izs:ize),INTENT(IN) :: Sapj, TColl_

    !! OUTPUT
    COMPLEX(dp), DIMENSION(ips:ipe,ijs:ije,ikys:ikye,ikxs:ikxe,izs:ize),INTENT(OUT) :: moments_rhs_

    INTEGER     :: p_int, j_int ! loops indices and polynom. degrees
    REAL(dp)    :: kx, ky, kperp2
    COMPLEX(dp) :: Tnapj, Tnapp2j, Tnapm2j, Tnapjp1, Tnapjm1 ! Terms from b x gradB and drives
    COMPLEX(dp) :: Tnapp1j, Tnapm1j, Tnapp1jm1, Tnapm1jm1 ! Terms from mirror force with non adiab moments_
    COMPLEX(dp) :: Tperp, Tpar, Tmir, Tphi, Tpsi
    COMPLEX(dp) :: Unapm1j, Unapm1jp1, Unapm1jm1 ! Terms from mirror force with adiab moments_
    COMPLEX(dp) :: i_kx,i_ky,phikykxz, psikykxz

     ! Measuring execution time
    CALL cpu_time(t0_rhs)

    ! Spatial loops
    zloop : DO  iz = izs,ize
      kxloop : DO ikx = ikxs,ikxe
        kx       = kxarray(ikx)   ! radial wavevector
        i_kx     = imagu * kx     ! radial derivative

        kyloop : DO iky = ikys,ikye
          ky     = kyarray(iky)   ! binormal wavevector
          i_ky   = imagu * ky     ! binormal derivative
          phikykxz = phi(iky,ikx,iz)! tmp phi value
          psikykxz = psi(iky,ikx,iz)! tmp psi value

          ! Kinetic loops
          jloop : DO ij = ijs, ije  ! This loop is from 1 to jmaxi+1
            j_int = jarray(ij)

            ploop : DO ip = ips, ipe  ! Hermite loop
              p_int = parray(ip)      ! Hermite degree
              eo    = MODULO(p_int,2) ! Indicates if we are on odd or even z grid
              kperp2= kparray(iky,ikx,iz,eo)**2

            IF((CLOS .NE. 1) .OR. (p_int+2*j_int .LE. dmaxe)) THEN
              !! Compute moments_ mixing terms
              Tperp = 0._dp; Tpar = 0._dp; Tmir = 0._dp
              ! Perpendicular dynamic
              ! term propto n^{p,j}
              Tnapj   = xnapj(ip,ij)* nadiab_moments(ip,ij,iky,ikx,iz)
              ! term propto n^{p+2,j}
              Tnapp2j = xnapp2j(ip) * nadiab_moments(ip+pp2,ij,iky,ikx,iz)
              ! term propto n^{p-2,j}
              Tnapm2j = xnapm2j(ip) * nadiab_moments(ip-pp2,ij,iky,ikx,iz)
              ! term propto n^{p,j+1}
              Tnapjp1 = xnapjp1(ij) * nadiab_moments(ip,ij+1,iky,ikx,iz)
              ! term propto n^{p,j-1}
              Tnapjm1 = xnapjm1(ij) * nadiab_moments(ip,ij-1,iky,ikx,iz)
              ! Tperp
              Tperp   = Tnapj + Tnapp2j + Tnapm2j + Tnapjp1 + Tnapjm1
              ! Parallel dynamic
              ! ddz derivative for Landau damping term
              Tpar      = xnapp1j(ip) * ddz_napj(ip+1,ij,iky,ikx,iz) &
                        + xnapm1j(ip) * ddz_napj(ip-1,ij,iky,ikx,iz)
              ! Mirror terms
              Tnapp1j   = ynapp1j  (ip,ij) * interp_napj(ip+1,ij  ,iky,ikx,iz)
              Tnapp1jm1 = ynapp1jm1(ip,ij) * interp_napj(ip+1,ij-1,iky,ikx,iz)
              Tnapm1j   = ynapm1j  (ip,ij) * interp_napj(ip-1,ij  ,iky,ikx,iz)
              Tnapm1jm1 = ynapm1jm1(ip,ij) * interp_napj(ip-1,ij-1,iky,ikx,iz)
              ! Trapping terms
              Unapm1j   = znapm1j  (ip,ij) * interp_napj(ip-1,ij  ,iky,ikx,iz)
              Unapm1jp1 = znapm1jp1(ip,ij) * interp_napj(ip-1,ij+1,iky,ikx,iz)
              Unapm1jm1 = znapm1jm1(ip,ij) * interp_napj(ip-1,ij-1,iky,ikx,iz)

              Tmir = Tnapp1j + Tnapp1jm1 + Tnapm1j + Tnapm1jm1 + Unapm1j + Unapm1jp1 + Unapm1jm1
              !! Electrical potential term
              IF ( p_int .LE. 2 ) THEN ! kronecker p0 p1 p2
                Tphi = (xphij  (ip,ij)*kernel(ij  ,iky,ikx,iz,eo) &
                      + xphijp1(ip,ij)*kernel(ij+1,iky,ikx,iz,eo) &
                      + xphijm1(ip,ij)*kernel(ij-1,iky,ikx,iz,eo))
              ELSE
                Tphi = 0._dp
              ENDIF

              !! Vector potential term
              IF ( (p_int .LE. 3) .AND. (p_int .GE. 1) ) THEN ! Kronecker p1 or p3
                Tpsi = (xpsij  (ip,ij)*kernel(ij  ,iky,ikx,iz,eo) &
                      + xpsijp1(ip,ij)*kernel(ij+1,iky,ikx,iz,eo) &
                      + xpsijm1(ip,ij)*kernel(ij-1,iky,ikx,iz,eo))
              ELSE
                Tpsi = 0._dp
              ENDIF

              !! Sum of all RHS terms
              moments_rhs_(ip,ij,iky,ikx,iz) = &
                  ! Perpendicular magnetic gradient/curvature effects
                  -imagu*Ckxky(iky,ikx,iz,eo) * Tperp&
                  ! Perpendicular pressure effects
                  -i_ky*beta*dpdx * Tperp&
                  ! Parallel coupling (Landau Damping)
                  -gradz_coeff(iz,eo) * Tpar &
                  ! Mirror term (parallel magnetic gradient)
                  -dBdz(iz,eo)*gradz_coeff(iz,eo) * Tmir&
                  ! Drives (density + temperature gradients)
                  -i_ky * (Tphi*phikykxz - Tpsi*psikykxz) &
                  ! Parallel drive term (should be negligible, test)
                  ! -Gamma_phipar(iz,eo)*Tphi*ddz_phi(iky,ikx,iz) &
                  ! Numerical Hermite hyperdiffusion (GX version)
                  -mu_p*diff_pe_coeff*p_int**4*moments_(ip,ij,iky,ikx,iz)&
                  ! Numerical Laguerre hyperdiffusion (GX version)
                  -mu_j*diff_je_coeff*j_int**4*moments_(ip,ij,iky,ikx,iz)&
                  ! Numerical perpendicular hyperdiffusion (totally artificial, for stability purpose)
                  -mu_x*diff_kx_coeff*kx**N_HD*moments_(ip,ij,iky,ikx,iz) &
                  -mu_y*diff_ky_coeff*ky**N_HD*moments_(ip,ij,iky,ikx,iz) &
                  ! Numerical parallel hyperdiffusion "mu_z*ddz**4"  see Pueschel 2010 (eq 25)
                  -mu_z*diff_dz_coeff*ddzND_napj(ip,ij,iky,ikx,iz) &
                  ! Collision term
                  +TColl_(ip,ij,iky,ikx,iz) &
                  ! Nonlinear term
                  -Gamma_NL(iz,eo)*Sapj(ip,ij,iky,ikx,iz)

                ! IF( (ip-4 .GT. 0) .AND. (num_procs_p .EQ. 1) ) &
                ! ! Numerical parallel velocity hyperdiffusion "+ dvpar4 g_a" see Pueschel 2010 (eq 33)
                ! ! (not used often so not parallelized)
                ! moments_rhs_(ip,ij,iky,ikx,iz) = &
                !   moments_rhs_(ip,ij,iky,ikx,iz) &
                !     + mu_p * moments_(ip-4,ij,iky,ikx,iz)

            ELSE
              moments_rhs_(ip,ij,iky,ikx,iz) = 0._dp
            ENDIF
            END DO ploop
          END DO jloop
        END DO kyloop
      END DO kxloop
    END DO zloop

    ! Execution time end
    CALL cpu_time(t1_rhs)
    tc_rhs = tc_rhs + (t1_rhs-t0_rhs)

  END SUBROUTINE moments_eq_rhs
  !_____________________________________________________________________________!
  !_____________________________________________________________________________!

END SUBROUTINE compute_moments_eq_rhs

SUBROUTINE add_Maxwellian_background_terms
  ! This routine is meant to add the terms rising from the magnetic operator,
  ! i.e. (B x GradB) Grad, applied on the background Maxwellian distribution
  ! (x_a + spar^2)(b x GradB) GradFaM
  ! It gives birth to kx=ky=0 sources terms (averages) that hit moments 00, 20,
  ! 40, 01,02, 21 with background gradient dependences.
  USE prec_const
  USE time_integration, ONLY : updatetlevel
  USE model,      ONLY: taue_qe, taui_qi, k_Ni, k_Ne, k_Ti, k_Te, KIN_E
  USE array,      ONLY: moments_rhs_e, moments_rhs_i
  USE grid,       ONLY: contains_kx0, contains_ky0, ikx_0, iky_0,&
                        ips_e,ipe_e,ijs_e,ije_e,ips_i,ipe_i,ijs_i,ije_i,&
                        zarray, izs,ize,&
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
                  +taue_qe * sinz(izs:ize) * (1.5_dp*k_Ne - 1.125_dp*k_Te)
            CASE(2) ! Na20 term
                moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                  +taue_qe * sinz(izs:ize) * (SQRT2*0.5_dp*k_Ne - 2.75_dp*k_Te)
            CASE(4) ! Na40 term
                moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                  +taue_qe * sinz(izs:ize) * SQRT6*0.75_dp*k_Te
            END SELECT
          CASE(1) ! j = 1
            SELECT CASE (ip-1)
            CASE(0) ! Na01 term
                moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                  -taue_qe * sinz(izs:ize) * (k_Ne + 3.5_dp*k_Te)
            CASE(2) ! Na21 term
                moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                  -taue_qe * sinz(izs:ize) * SQRT2*k_Te
            END SELECT
          CASE(2) ! j = 2
            SELECT CASE (ip-1)
            CASE(0) ! Na02 term
                moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_e(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                  +taue_qe * sinz(izs:ize) * 2._dp*k_Te
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
                +taui_qi * sinz(izs:ize) * (1.5_dp*k_Ni - 1.125_dp*k_Ti)
          CASE(2) ! Na20 term
              moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                +taui_qi * sinz(izs:ize) * (SQRT2*0.5_dp*k_Ni - 2.75_dp*k_Ti)
          CASE(4) ! Na40 term
              moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                +taui_qi * sinz(izs:ize) * SQRT6*0.75_dp*k_Ti
          END SELECT
        CASE(1) ! j = 1
          SELECT CASE (ip-1)
          CASE(0) ! Na01 term
              moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                -taui_qi * sinz(izs:ize) * (k_Ni + 3.5_dp*k_Ti)
          CASE(2) ! Na21 term
              moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                -taui_qi * sinz(izs:ize) * SQRT2*k_Ti
          END SELECT
        CASE(2) ! j = 2
          SELECT CASE (ip-1)
          CASE(0) ! Na02 term
              moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel) = moments_rhs_i(ip,ij,iky_0,ikx_0,izs:ize,updatetlevel)&
                +taui_qi * sinz(izs:ize) * 2._dp*k_Ti
          END SELECT
        END SELECT
      ENDDO
    ENDDO

  ENDIF

END SUBROUTINE

END MODULE moments_eq_rhs
