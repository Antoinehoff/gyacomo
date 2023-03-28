MODULE moments_eq_rhs
  IMPLICIT NONE
  PUBLIC :: compute_moments_eq_rhs
CONTAINS

SUBROUTINE compute_moments_eq_rhs
  USE model
  USE array
  USE fields
  USE grid,       ONLY: local_na, local_np, local_nj, local_nkx, local_nky, local_nz,&
                        nzgrid,pp2,ngp,ngj,ngz,dmax,&
                        diff_dz_coeff,diff_kx_coeff,diff_ky_coeff,diff_p_coeff,diff_j_coeff,&
                        parray,jarray,kxarray, kyarray, kparray
  USE basic
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
              IF((CLOS .NE. 1) .OR. (p_int +2*j_int .LE. dmax)) THEN ! for the closure scheme
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
                  Dphi =i_ky*( xphij  (ia,ip,ij)*kernel(ia,iji  ,iky,ikx,izi,eo) &
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

  !     print*,'sumabs moments i', SUM(ABS(moments(1,:,:,:,:,:,updatetlevel)))
  !     print*,'moment rhs i 221', moments_rhs(1,1,1,2,2,1,updatetlevel)
  !     ! print*,'sum real nadiabe', SUM(REAL(nadiab_moments(2,:,:,:,:,:)))
  !     print*,'sumreal momrhs i', SUM(REAL(moments_rhs(1,:,:,:,:,:,:)))
  !     print*,'sumimag momrhs i', SUM(IMAG(moments_rhs(1,:,:,:,:,:,:)))
  !     print*,'sumreal phi     ', SUM(REAL(phi(:,:,(1+ngz/2):(local_nz+ngz/2))))
  !     print*,'sumimag phi     ', SUM(IMAG(phi(:,:,(1+ngz/2):(local_nz+ngz/2))))
  !     print*,'phi(2,2,1)      ', phi(2,2,1+ngz/2)
  !     print*,'sumreal ddznipj ', SUM(REAL(ddz_napj(1,:,:,:,:,:)))
  !     print*,'sumimag ddznipj ', SUM(IMAG(ddz_napj(1,:,:,:,:,:)))
  !     print*,'        ddznipj  ',(ddz_napj(1,2+ngp/2,2+ngj/2,2,2,1))
  !     print*,'      ddzNDnipj  ',SUM(REAL(ddzND_Napj(1,:,:,:,:,:)))
  !     print*,'sumreal Capj    ', SUM(REAL(Capj(1,:,:,:,:,:)))
  !     print*,'---'
  !     IF(updatetlevel .EQ. 4) stop
  ! stop

END SUBROUTINE compute_moments_eq_rhs

! SUBROUTINE add_Maxwellian_background_terms
!   ! This routine is meant to add the terms rising from the magnetic operator,
!   ! i.e. (B x k_gB) Grad, applied on the background Maxwellian distribution
!   ! (x_a + spar^2)(b x k_gB) GradFaM
!   ! It gives birth to kx=ky=0 sources terms (averages) that hit moments_ 00, 20,
!   ! 40, 01,02, 21 with background gradient dependences.
!   USE prec_const
!   USE time_integration, ONLY : updatetlevel
!   USE species,          ONLY: tau_q, k_N, k_T
!   USE array,            ONLY: moments_rhs
!   USE grid,             ONLY: contains_kx0, contains_ky0, ikx0, iky0,&
!                               ia,ias,iae,ip,ips,ipe, ij,ijs,ije, zarray,izs,ize
!   IMPLICIT NONE
!   real(xp), DIMENSION(izs:ize) :: sinz
!
!   sinz(izs:ize) = SIN(zarray(izs:ize,0))
!
!   IF(contains_kx0 .AND. contains_ky0) THEN
!     DO ia = ias, iae
!       DO ip = ips,ipe
!         DO ij = ijs,ije
!           SELECT CASE(ij-1)
!           CASE(0) ! j = 0
!             SELECT CASE (ip-1)
!             CASE(0) ! Na00 term
!                 moments_rhs(ia,ip,ij,iky0,ikx0,izs:ize,updatetlevel) = moments_rhs(ia,ip,ij,iky0,ikx0,izs:ize,updatetlevel)&
!                   +tau_q(ia) * sinz(izs:ize) * (1.5_xp*k_N(ia) - 1.125_xp*k_T(ia))
!             CASE(2) ! Na20 term
!                 moments_rhs(ia,ip,ij,iky0,ikx0,izs:ize,updatetlevel) = moments_rhs(ia,ip,ij,iky0,ikx0,izs:ize,updatetlevel)&
!                   +tau_q(ia) * sinz(izs:ize) * (SQRT2*0.5_xp*k_N(ia) - 2.75_xp*k_T(ia))
!             CASE(4) ! Na40 term
!                 moments_rhs(ia,ip,ij,iky0,ikx0,izs:ize,updatetlevel) = moments_rhs(ia,ip,ij,iky0,ikx0,izs:ize,updatetlevel)&
!                   +tau_q(ia) * sinz(izs:ize) * SQRT6*0.75_xp*k_T(ia)
!             END SELECT
!           CASE(1) ! j = 1
!             SELECT CASE (ip-1)
!             CASE(0) ! Na01 term
!                 moments_rhs(ia,ip,ij,iky0,ikx0,izs:ize,updatetlevel) = moments_rhs(ia,ip,ij,iky0,ikx0,izs:ize,updatetlevel)&
!                   -tau_q(ia) * sinz(izs:ize) * (k_N(ia) + 3.5_xp*k_T(ia))
!             CASE(2) ! Na21 term
!                 moments_rhs(ia,ip,ij,iky0,ikx0,izs:ize,updatetlevel) = moments_rhs(ia,ip,ij,iky0,ikx0,izs:ize,updatetlevel)&
!                   -tau_q(ia) * sinz(izs:ize) * SQRT2*k_T(ia)
!             END SELECT
!           CASE(2) ! j = 2
!             SELECT CASE (ip-1)
!             CASE(0) ! Na02 term
!                 moments_rhs(ia,ip,ij,iky0,ikx0,izs:ize,updatetlevel) = moments_rhs(ia,ip,ij,iky0,ikx0,izs:ize,updatetlevel)&
!                   +tau_q(ia) * sinz(izs:ize) * 2._xp*k_T(ia)
!             END SELECT
!           END SELECT
!         ENDDO
!       ENDDO
!     ENDDO
!   ENDIF
!
! END SUBROUTINE

END MODULE moments_eq_rhs
