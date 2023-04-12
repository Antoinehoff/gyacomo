MODULE moments_eq_rhs
  IMPLICIT NONE
  PUBLIC :: compute_moments_eq_rhs
CONTAINS

  SUBROUTINE compute_moments_eq_rhs
    USE model
    USE array
    USE fields
    USE grid,       ONLY: local_na, local_np, local_nj, local_nkx, local_nky, local_nz,&
                          nzgrid,pp2,ngp,ngj,ngz,&
                          diff_dz_coeff,diff_kx_coeff,diff_ky_coeff,diff_p_coeff,diff_j_coeff,&
                          parray,jarray,kxarray, kyarray, kparray
    USE basic
    USE closure,    ONLY: evolve_mom
    USE prec_const
    USE collision
    USE time_integration
    ! USE species, ONLY: xpdx
    USE geometry, ONLY: gradz_coeff, dlnBdz, Ckxky!, Gamma_phipar
    USE calculus, ONLY: interp_z, grad_z, grad_z2
    ! USE species,  ONLY: xpdx
    IMPLICIT NONE
    INTEGER     :: ia, iz, iky,  ikx, ip ,ij, eo ! counters
    INTEGER     :: izi,ipi,iji ! interior points counters
    INTEGER     :: p_int, j_int ! loops indices and polynom. degrees
    REAL(xp)    :: kx, ky, kperp2
    COMPLEX(xp) :: Tnapj, Tnapp2j, Tnapm2j, Tnapjp1, Tnapjm1 ! Terms from b x gradB and drives
    COMPLEX(xp) :: Tnapp1j, Tnapm1j, Tnapp1jm1, Tnapm1jm1 ! Terms from mirror force with non adiab moments_
    COMPLEX(xp) :: Ldamp, Fmir
    COMPLEX(xp) :: Mperp, Mpara, Dphi, Dpsi
    COMPLEX(xp) :: Unapm1j, Unapm1jp1, Unapm1jm1 ! Terms from mirror force with adiab moments_
    COMPLEX(xp) :: i_kx,i_ky
    COMPLEX(xp) :: Napj, RHS, phikykxz, psikykxz
    ! Spatial loops
    z:DO  iz = 1,local_nz
      izi = iz + ngz/2
      x:DO ikx = 1,local_nkx
        kx       = kxarray(ikx)                     ! radial wavevector
        i_kx     = imagu * kx                       ! radial derivative
        y:DO iky = 1,local_nky
          ky       = kyarray(iky)                     ! binormal wavevector
          i_ky     = imagu * ky                       ! binormal derivative
          phikykxz = phi(iky,ikx,izi)
          psikykxz = psi(iky,ikx,izi)
          ! Kinetic loops
          j:DO ij = 1,local_nj               ! This loop is from 1 to jmaxi+1
            iji   = ij+ngj/2
            j_int = jarray(iji)
            p:DO ip = 1,local_np             ! Hermite loop
              ipi   = ip+ngp/2
              p_int = parray(ipi)                   ! Hermite degree
              eo    = min(nzgrid,MODULO(p_int,2)+1) ! Indicates if we are on odd or even z grid
              kperp2= kparray(iky,ikx,izi,eo)**2
              ! Species loop
              a:DO ia = 1,local_na
                Napj = moments(ia,ipi,iji,iky,ikx,izi,updatetlevel)
                RHS = 0._xp
                IF(evolve_mom(ipi,iji)) THEN ! for the closure scheme
                  !! Compute moments_ mixing terms
                  ! Perpendicular dynamic
                  ! term propto n^{p,j}
                  Tnapj   = xnapj(ia,ip,ij)* nadiab_moments(ia,ipi,    iji,  iky,ikx,izi)
                  ! term propto n^{p+2,j}
                  Tnapp2j = xnapp2j(ia,ip) * nadiab_moments(ia,ipi+pp2,iji,  iky,ikx,izi)
                  ! term propto n^{p-2,j}
                  Tnapm2j = xnapm2j(ia,ip) * nadiab_moments(ia,ipi-pp2,iji,  iky,ikx,izi)
                  ! term propto n^{p,j+1}
                  Tnapjp1 = xnapjp1(ia,ij) * nadiab_moments(ia,ipi,    iji+1,iky,ikx,izi)
                  ! term propto n^{p,j-1}
                  Tnapjm1 = xnapjm1(ia,ij) * nadiab_moments(ia,ipi,    iji-1,iky,ikx,izi)
                  ! Perpendicular magnetic term (curvature and gradient drifts)
                  Mperp   = imagu*Ckxky(iky,ikx,izi,eo)&
                            *(Tnapj + Tnapp2j + Tnapm2j + Tnapjp1 + Tnapjm1)
                  ! Parallel dynamic
                  ! ddz derivative for Landau damping term
                  Ldamp     = xnapp1j(ia,ip) * ddz_napj(ia,ipi+1,iji,iky,ikx,iz) &
                            + xnapm1j(ia,ip) * ddz_napj(ia,ipi-1,iji,iky,ikx,iz)
                  ! Mirror terms
                  Tnapp1j   = ynapp1j  (ia,ip,ij) * interp_napj(ia,ipi+1,iji  ,iky,ikx,iz)
                  Tnapp1jm1 = ynapp1jm1(ia,ip,ij) * interp_napj(ia,ipi+1,iji-1,iky,ikx,iz)
                  Tnapm1j   = ynapm1j  (ia,ip,ij) * interp_napj(ia,ipi-1,iji  ,iky,ikx,iz)
                  Tnapm1jm1 = ynapm1jm1(ia,ip,ij) * interp_napj(ia,ipi-1,iji-1,iky,ikx,iz)
                  ! Trapping terms
                  Unapm1j   = znapm1j  (ia,ip,ij) * interp_napj(ia,ipi-1,iji  ,iky,ikx,iz)
                  Unapm1jp1 = znapm1jp1(ia,ip,ij) * interp_napj(ia,ipi-1,iji+1,iky,ikx,iz)
                  Unapm1jm1 = znapm1jm1(ia,ip,ij) * interp_napj(ia,ipi-1,iji-1,iky,ikx,iz)
                  ! sum the parallel forces
                  Fmir = dlnBdz(izi,eo)*(Tnapp1j + Tnapp1jm1 + Tnapm1j + Tnapm1jm1 +&
                                        Unapm1j + Unapm1jp1 + Unapm1jm1)
                  ! Parallel magnetic term (Landau damping and the mirror force)
                  Mpara = gradz_coeff(izi,eo)*(Ldamp + Fmir)
                  !! Electrical potential term
                  IF ( p_int  .LE. 2 ) THEN ! kronecker p0 p1 p2
                    Dphi =i_ky*( xphij(ia,ip,ij)*kernel(ia,iji  ,iky,ikx,izi,eo) &
                              +xphijp1(ia,ip,ij)*kernel(ia,iji+1,iky,ikx,izi,eo) &
                              +xphijm1(ia,ip,ij)*kernel(ia,iji-1,iky,ikx,izi,eo) )*phikykxz
                  ELSE
                    Dphi = 0._xp
                  ENDIF
                  !! Vector potential term
                  IF ( (p_int  .LE. 3) .AND. (p_int  .GE. 1) ) THEN ! Kronecker p1 or p3
                    Dpsi =-i_ky*( xpsij  (ia,ip,ij)*kernel(ia,iji  ,iky,ikx,izi,eo) &
                                +xpsijp1(ia,ip,ij)*kernel(ia,iji+1,iky,ikx,izi,eo) &
                                +xpsijm1(ia,ip,ij)*kernel(ia,iji-1,iky,ikx,izi,eo))*psikykxz
                  ELSE
                    Dpsi = 0._xp
                  ENDIF
                  !! Sum of all RHS terms
                  RHS = &
                      ! Nonlinear term Sapj_ = {phi,f}
                      - Sapj(ia,ip,ij,iky,ikx,iz) &
                      ! Perpendicular magnetic term
                      - Mperp &
                      ! Parallel magnetic term
                      - Mpara &
                      ! Drives (density + temperature gradients)
                      - (Dphi + Dpsi) &
                      ! Collision term
                      + Capj(ia,ip,ij,iky,ikx,iz) &
                      ! Perpendicular pressure effects (electromagnetic term) (TO CHECK)
                      ! - i_ky*beta*xpdx(ia) * (Tnapj + Tnapp2j + Tnapm2j + Tnapjp1 + Tnapjm1)&
                      ! Parallel drive term (should be negligible, to test)
                      ! !! -Gamma_phipar(iz,eo)*Tphi*ddz_phi(iky,ikx,iz) &
                      ! Numerical perpendicular hyperdiffusion
                      -mu_x*diff_kx_coeff*kx**N_HD*Napj &
                      -mu_y*diff_ky_coeff*ky**N_HD*Napj &
                      ! ! ! Numerical parallel hyperdiffusion "mu_z*ddz**4"  see Pueschel 2010 (eq 25)
                      -mu_z*diff_dz_coeff*ddzND_napj(ia,ipi,iji,iky,ikx,iz)
                  !! Velocity space dissipation (should be implemented somewhere else)
                  SELECT CASE(HYP_V)
                  CASE('hypcoll') ! GX like Hermite hypercollisions see Mandell et al. 2023 (eq 3.23), unadvised to use it
                    IF (p_int  .GT. 2)  &
                      RHS = RHS - mu_p*diff_p_coeff*p_int**6*Napj
                    IF (j_int .GT. 1)  &
                      RHS = RHS - mu_j*diff_j_coeff*j_int**6*Napj
                  CASE('dvpar4')
                    ! fourth order numerical diffusion in vpar
                    IF(p_int  .GE. 4) &
                    ! Numerical parallel velocity hyperdiffusion "+ dvpar4 g_a" see Pueschel 2010 (eq 33)
                    ! (not used often so not parallelized)
                    RHS = RHS + mu_p*dv4_Hp_coeff(p_int)*moments(ia,ipi-4,iji,iky,ikx,izi,updatetlevel)
                    ! + dummy Laguerre diff
                    IF (j_int .GE. 2)  &
                      RHS = RHS - mu_j*diff_j_coeff*j_int**6*Napj
                  CASE DEFAULT
                  END SELECT
                ELSE
                  RHS = 0._xp
                ENDIF
                !! Put RHS in the array
                moments_rhs(ia,ip,ij,iky,ikx,iz,updatetlevel) = RHS
              END DO a
            END DO p
          END DO j
        END DO y
      END DO x
    END DO z
  END SUBROUTINE compute_moments_eq_rhs

END MODULE moments_eq_rhs
