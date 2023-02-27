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
  USE geometry, ONLY: gradz_coeff, dlnBdz, Ckxky!, Gamma_phipar
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
                     TColl_i, ddzND_Nipj, diff_pi_coeff, diff_ji_coeff,&
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
                     TColl_e, ddzND_Nepj, diff_pe_coeff, diff_je_coeff,&
                     moments_rhs_e(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,ikxs:ikxe,izs:ize,updatetlevel))

  CONTAINS
  !_____________________________________________________________________________!
  !_____________________________________________________________________________!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! moments_ RHS computation !!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This routine assemble the RHS of the moment hierarchy equations. It uses
  ! linear coefficients that are stored in arrays (xn*, yn* and zn*) computed in
  ! numerics_mod.F90. Otherwise it simply adds the collision term TColl_ that is
  ! computed in collision_mod.F90 and the nonlinear term Sapj_ computed in
  ! nonlinear_mod.F90.
  ! All arguments of the subroutines are inputs only except the last one,
  ! moments_rhs_ that will contain the sum of every terms in the RHS.
  !_____________________________________________________________________________!
  SUBROUTINE moments_eq_rhs(ips_,ipe_,ipgs_,ipge_,ijs_,ije_,ijgs_,ijge_,jarray_,parray_,&
                   xnapj_, xnapp2j_, xnapm2j_, xnapjp1_, xnapjm1_, xnapp1j_, xnapm1j_,&
                   ynapp1j_, ynapp1jm1_, ynapm1j_, ynapm1jm1_, &
                   znapm1j_, znapm1jp1_, znapm1jm1_, &
                   xphij_, xphijp1_, xphijm1_, xpsij_, xpsijp1_, xpsijm1_,&
                   kernel_, nadiab_moments_, ddz_napj_, interp_napj_, Sapj_,&
                   moments_, TColl_, ddzND_napj_, diff_p_coeff_, diff_j_coeff_, moments_rhs_)

    IMPLICIT NONE
    !! INPUTS
    INTEGER, INTENT(IN) :: ips_, ipe_, ipgs_, ipge_, ijs_, ije_, ijgs_, ijge_
    INTEGER,  DIMENSION(ips_:ipe_), INTENT(IN) :: parray_
    INTEGER,  DIMENSION(ijs_:ije_), INTENT(IN) :: jarray_
    REAL(dp), DIMENSION(ips_:ipe_,ijs_:ije_), INTENT(IN) :: xnapj_
    REAL(dp), DIMENSION(ips_:ipe_),         INTENT(IN) :: xnapp2j_, xnapm2j_
    REAL(dp), DIMENSION(ijs_:ije_),         INTENT(IN) :: xnapjp1_, xnapjm1_
    REAL(dp), DIMENSION(ips_:ipe_),         INTENT(IN) :: xnapp1j_, xnapm1j_
    REAL(dp), DIMENSION(ips_:ipe_,ijs_:ije_), INTENT(IN) :: ynapp1j_, ynapp1jm1_, ynapm1j_, ynapm1jm1_
    REAL(dp), DIMENSION(ips_:ipe_,ijs_:ije_), INTENT(IN) :: znapm1j_, znapm1jp1_, znapm1jm1_
    REAL(dp), DIMENSION(ips_:ipe_,ijs_:ije_), INTENT(IN) :: xphij_, xphijp1_, xphijm1_
    REAL(dp), DIMENSION(ips_:ipe_,ijs_:ije_), INTENT(IN) :: xpsij_, xpsijp1_, xpsijm1_

    REAL(dp), DIMENSION(ijgs_:ijge_,ikys:ikye,ikxs:ikxe,izgs:izge,0:1),INTENT(IN) :: kernel_

    COMPLEX(dp), DIMENSION(ipgs_:ipge_,ijgs_:ijge_,ikys:ikye,ikxs:ikxe,izgs:izge),INTENT(IN) :: nadiab_moments_
    COMPLEX(dp), DIMENSION(ipgs_:ipge_,ijgs_:ijge_,ikys:ikye,ikxs:ikxe,izgs:izge),INTENT(IN) :: ddz_napj_
    COMPLEX(dp), DIMENSION(ipgs_:ipge_,ijgs_:ijge_,ikys:ikye,ikxs:ikxe,izgs:izge),INTENT(IN) :: interp_napj_
    COMPLEX(dp), DIMENSION( ips_:ipe_,  ijs_:ije_, ikys:ikye,ikxs:ikxe, izs:ize), INTENT(IN) :: Sapj_
    COMPLEX(dp), DIMENSION(ipgs_:ipge_,ijgs_:ijge_,ikys:ikye,ikxs:ikxe,izgs:izge),INTENT(IN) :: moments_
    COMPLEX(dp), DIMENSION( ips_:ipe_,  ijs_:ije_, ikys:ikye,ikxs:ikxe, izs:ize), INTENT(IN) :: TColl_
    COMPLEX(dp), DIMENSION(ipgs_:ipge_,ijgs_:ijge_,ikys:ikye,ikxs:ikxe,izgs:izge),INTENT(IN) :: ddzND_napj_
    REAL(dp),    INTENT(IN) :: diff_p_coeff_, diff_j_coeff_
    !! OUTPUT
    COMPLEX(dp), DIMENSION( ips_:ipe_,  ijs_:ije_, ikys:ikye,ikxs:ikxe, izs:ize),INTENT(OUT) :: moments_rhs_

    INTEGER     :: p_int, j_int ! loops indices and polynom. degrees
    REAL(dp)    :: kx, ky, kperp2
    COMPLEX(dp) :: Tnapj, Tnapp2j, Tnapm2j, Tnapjp1, Tnapjm1 ! Terms from b x gradB and drives
    COMPLEX(dp) :: Tnapp1j, Tnapm1j, Tnapp1jm1, Tnapm1jm1 ! Terms from mirror force with non adiab moments_
    COMPLEX(dp) :: Tpar, Tmir, Tphi, Tpsi
    COMPLEX(dp) :: Mperp, Mpara, Dphi, Dpsi
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
          psikykxz = psi(iky,ikx,iz)! tmp psi value

          ! Kinetic loops
          jloop : DO ij = ijs_, ije_  ! This loop is from 1 to jmaxi+1
            j_int = jarray_(ij)

            ploop : DO ip = ips_, ipe_  ! Hermite loop
              p_int = parray_(ip)      ! Hermite degree
              eo    = MODULO(p_int,2) ! Indicates if we are on odd or even z grid
              kperp2= kparray(iky,ikx,iz,eo)**2

            IF((CLOS .NE. 1) .OR. (p_int+2*j_int .LE. dmaxe)) THEN ! for the closure scheme
              !! Compute moments_ mixing terms
              ! Perpendicular dynamic
              ! term propto n^{p,j}
              Tnapj   = xnapj_(ip,ij)* nadiab_moments_(ip,ij,iky,ikx,iz)
              ! term propto n^{p+2,j}
              Tnapp2j = xnapp2j_(ip) * nadiab_moments_(ip+pp2,ij,iky,ikx,iz)
              ! term propto n^{p-2,j}
              Tnapm2j = xnapm2j_(ip) * nadiab_moments_(ip-pp2,ij,iky,ikx,iz)
              ! term propto n^{p,j+1}
              Tnapjp1 = xnapjp1_(ij) * nadiab_moments_(ip,ij+1,iky,ikx,iz)
              ! term propto n^{p,j-1}
              Tnapjm1 = xnapjm1_(ij) * nadiab_moments_(ip,ij-1,iky,ikx,iz)
              ! Perpendicular magnetic term (curvature and gradient drifts)
              Mperp   = imagu*Ckxky(iky,ikx,iz,eo)*(Tnapj + Tnapp2j + Tnapm2j + Tnapjp1 + Tnapjm1)

              ! Parallel dynamic
              ! ddz derivative for Landau damping term
              Tpar      = xnapp1j_(ip) * ddz_napj_(ip+1,ij,iky,ikx,iz) &
                        + xnapm1j_(ip) * ddz_napj_(ip-1,ij,iky,ikx,iz)
              ! Mirror terms
              Tnapp1j   = ynapp1j_  (ip,ij) * interp_napj_(ip+1,ij  ,iky,ikx,iz)
              Tnapp1jm1 = ynapp1jm1_(ip,ij) * interp_napj_(ip+1,ij-1,iky,ikx,iz)
              Tnapm1j   = ynapm1j_  (ip,ij) * interp_napj_(ip-1,ij  ,iky,ikx,iz)
              Tnapm1jm1 = ynapm1jm1_(ip,ij) * interp_napj_(ip-1,ij-1,iky,ikx,iz)
              ! Trapping terms
              Unapm1j   = znapm1j_  (ip,ij) * interp_napj_(ip-1,ij  ,iky,ikx,iz)
              Unapm1jp1 = znapm1jp1_(ip,ij) * interp_napj_(ip-1,ij+1,iky,ikx,iz)
              Unapm1jm1 = znapm1jm1_(ip,ij) * interp_napj_(ip-1,ij-1,iky,ikx,iz)

              Tmir = dlnBdz(iz,eo)*(Tnapp1j + Tnapp1jm1 + Tnapm1j + Tnapm1jm1 +&
                                    Unapm1j + Unapm1jp1 + Unapm1jm1)
              ! Parallel magnetic term (Landau damping and the mirror force)
              Mpara = gradz_coeff(iz,eo)*(Tpar + Tmir)
              !! Electrical potential term
              IF ( p_int .LE. 2 ) THEN ! kronecker p0 p1 p2
                Dphi =i_ky*( xphij_  (ip,ij)*kernel_(ij  ,iky,ikx,iz,eo) &
                            +xphijp1_(ip,ij)*kernel_(ij+1,iky,ikx,iz,eo) &
                            +xphijm1_(ip,ij)*kernel_(ij-1,iky,ikx,iz,eo) )*phi(iky,ikx,iz)
              ELSE
                Tphi = 0._dp
              ENDIF

              !! Vector potential term
              IF ( (p_int .LE. 3) .AND. (p_int .GE. 1) ) THEN ! Kronecker p1 or p3
                Dpsi =-i_ky*( xpsij_  (ip,ij)*kernel_(ij  ,iky,ikx,iz,eo) &
                             +xpsijp1_(ip,ij)*kernel_(ij+1,iky,ikx,iz,eo) &
                             +xpsijm1_(ip,ij)*kernel_(ij-1,iky,ikx,iz,eo))*psi(iky,ikx,iz)
              ELSE
                Dpsi = 0._dp
              ENDIF

              !! Sum of all RHS terms
              moments_rhs_(ip,ij,iky,ikx,iz) = &
                  ! Nonlinear term Sapj_ = {phi,f}
                  - Sapj_(ip,ij,iky,ikx,iz) &
                  ! Perpendicular magnetic term
                  - Mperp &
                  ! Parallel magnetic term
                  - Mpara &
                  ! Drives (density + temperature gradients)
                  - (Dphi + Dpsi) &
                  ! Collision term
                  + TColl_(ip,ij,iky,ikx,iz) &
                  ! Perpendicular pressure effects (electromagnetic term) (TO CHECK)
                  - i_ky*beta*dpdx * (Tnapj + Tnapp2j + Tnapm2j + Tnapjp1 + Tnapjm1)&
                  ! Parallel drive term (should be negligible, to test)
                  ! -Gamma_phipar(iz,eo)*Tphi*ddz_phi(iky,ikx,iz) &
                  ! Numerical perpendicular hyperdiffusion
                  -mu_x*diff_kx_coeff*kx**N_HD*moments_(ip,ij,iky,ikx,iz) &
                  -mu_y*diff_ky_coeff*ky**N_HD*moments_(ip,ij,iky,ikx,iz) &
                  ! Numerical parallel hyperdiffusion "mu_z*ddz**4"  see Pueschel 2010 (eq 25)
                  -mu_z*diff_dz_coeff*ddzND_napj_(ip,ij,iky,ikx,iz)

                  !! Velocity space dissipation (should be implemented somewhere else)
                  SELECT CASE(HYP_V)
                  CASE('hypcoll') ! GX like Hermite hypercollisions see Mandell et al. 2023 (eq 3.23), unadvised to use it
                    IF (p_int .GT. 2)  &
                      moments_rhs_(ip,ij,iky,ikx,iz) = &
                        moments_rhs_(ip,ij,iky,ikx,iz) - mu_p*diff_pe_coeff*p_int**6*moments_(ip,ij,iky,ikx,iz)
                    IF (j_int .GT. 1)  &
                      moments_rhs_(ip,ij,iky,ikx,iz) = &
                        moments_rhs_(ip,ij,iky,ikx,iz) - mu_j*diff_je_coeff*j_int**6*moments_(ip,ij,iky,ikx,iz)
                  CASE('dvpar4')
                    ! fourth order numerical diffusion in vpar
                    IF(ip-4 .GT. 0) &
                    ! Numerical parallel velocity hyperdiffusion "+ dvpar4 g_a" see Pueschel 2010 (eq 33)
                    ! (not used often so not parallelized)
                    moments_rhs_(ip,ij,iky,ikx,iz) = &
                      moments_rhs_(ip,ij,iky,ikx,iz) &
                        + mu_p*dv4_Hp_coeff(p_int)*moments_(ip-4,ij,iky,ikx,iz)
                    ! + dummy Laguerre diff
                    IF (j_int .GT. 1)  &
                      moments_rhs_(ip,ij,iky,ikx,iz) = &
                        moments_rhs_(ip,ij,iky,ikx,iz) - mu_j*diff_je_coeff*j_int**6*moments_(ip,ij,iky,ikx,iz)
                  CASE DEFAULT
                  END SELECT
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
  ! i.e. (B x k_gB) Grad, applied on the background Maxwellian distribution
  ! (x_a + spar^2)(b x k_gB) GradFaM
  ! It gives birth to kx=ky=0 sources terms (averages) that hit moments_ 00, 20,
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
