!! This source has been adapted from the GENE code https://genecode.org/ !!
!>Provides a circular concentric geometry
MODULE circular
  USE prec_const
  USE basic
  USE parallel, ONLY: my_id, num_procs_z, nbr_U, nbr_D, comm0
  ! use coordinates,only: gcoor, get_dzprimedz
  USE grid, ONLY: local_Nky, local_Nkx, local_nz, ngz, nzgrid, kyarray, kxarray, zarray, total_nz, local_nz_offset, iodd, ieven
  ! use discretization
  USE lagrange_interpolation

  implicit none
  public:: get_circ
  private

! SUBROUTINES DEFINITIONS
CONTAINS

  SUBROUTINE get_circ(trpeps,major_R,major_Z,q0,shat,&
    C_y,C_xy,Cyq0_x0,gxx_,gxy_,gxz_,gyy_,gyz_,gzz_,dBdx_,dBdy_,dBdz_,&
    Bfield_,jacobian_,R_hat_,Z_hat_,dxdR_,dxdZ_)
    !!!!!!!!!!!!!!!! GYACOMO INTERFACE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(xp), INTENT(INOUT) :: trpeps          ! eps in gyacomo (inverse aspect ratio)
    real(xp), INTENT(INOUT) :: major_R         ! major radius
    real(xp), INTENT(INOUT) :: major_Z         ! major Z
    real(xp), INTENT(INOUT) :: q0              ! safetyfactor
    real(xp), INTENT(INOUT) :: shat            ! safetyfactor
    real(xp), INTENT(INOUT) :: C_y, C_xy, Cyq0_x0
    real(xp), dimension(local_nz+ngz,nzgrid), INTENT(OUT) :: &
                                              gxx_,gxy_,gxz_,gyy_,gyz_,gzz_,&
                                              dBdx_,dBdy_,dBdz_,Bfield_,jacobian_,&
                                              R_hat_,Z_hat_,dxdR_,dxdZ_
    INTEGER :: iz,eo
    ! No parameter in gyacomo yet
    real(xp) :: SIGN_Ip=1._xp ! current sign (only normal current)
    !!!!!!!!!!!!!! END GYACOMO INTERFACE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    REAL(xp)    :: g11,g12,g22,g33 ! normalized metric in (PSI,CHI,PHI)
    REAL(xp)    :: B,dBdr_c,dBdchi ! J_PCP : J_PSI,CHI,P
    ! intermediate variables, t stands for theta the poloidal angle.
    REAL(xp)    :: qbar,dpsidr,dchidr,dchidt,dBdr,dBdt,cost,sint
    REAL(xp)    :: fac2, dCydr_Cy
    REAL(xp)    :: dxdr, dqdr, dqdx
    REAL(xp)    :: dqbardr, minor_r, r0, z

    minor_r  = major_R
    r0       = major_R*trpeps
    qbar     = abs(q0)*sqrt(1-trpeps**2)
    C_y      = abs(r0/q0)*SIGN_Ip
    dxdr     = 1._xp                ! x=r
    dpsidr   = r0/qbar                ! 1/(Bref*Lref)*dpsi/dr
    ! dqdx     = shat*q0/minor_r/dxdr ! AH: added
    dqdx     = shat/C_y  ! AH: added
    dqdr     = dqdx*dxdr/minor_r
    dqbardr  = sqrt(1._xp-trpeps**2)*(abs(dqdr)-abs(q0)*trpeps/major_R/(1-trpeps**2))
    dCydr_Cy = 0._xp
    C_xy     = dpsidr/(dxdr*abs(C_y))
    Cyq0_x0  = C_y*q0/r0
    fac2 = (dCydr_Cy+dqdr/q0)

    DO eo = 1,nzgrid
      do iz=1,local_nz+ngz
        z = zarray(iz,eo) 
        ! z dep. intermediate variables
        cost = (cos(z) - trpeps)/(1._xp-trpeps*cos(z))
        sint = sqrt(1._xp-trpeps**2)*sin(SIGN_Ip*z)/(1._xp-trpeps*cos(z))
        ! derivatives with respect to r,theta
        dchidr =  -1._xp/major_R*sin(z)/(1._xp-trpeps**2)
        dchidt = SIGN_Ip * sqrt(1-trpeps**2)/(1+trpeps*cost)

        B    = 1._xp/(1._xp+trpeps*cost)*sqrt(1.+(trpeps/qbar)**2)
        dBdr = B/major_R*(-cost/(1._xp+trpeps*cost) + &
               + trpeps/(trpeps**2+qbar**2)*(1._xp - major_R*trpeps*dqbardr/qbar))
        dBdt = trpeps*sint*B / (1._xp+trpeps*cost)

        ! metric in (r,CHI,PHI)
        g11 = 1._xp
        g22 = dchidr**2+1._xp/(major_R*trpeps)**2*dchidt**2
        g12 = dchidr
        g33 = 1._xp/major_R**2*1._xp/(1._xp+trpeps*cost)**2

        ! magnetic field derivatives in
        dBdchi = dBdt/dchidt
        dBdr_c = 1._xp/g11*(dBdr-g12*dBdt/dchidt)

        ! metric in (x,y,z)
        gxx_(iz,eo) = (dxdr)**2*g11
        gyy_(iz,eo) = (C_y*q0)**2 * ((fac2*z)**2*g11 +2._xp*fac2*z*g12+g22) + C_y**2*g33
        gxy_(iz,eo) = dxdr*C_y*SIGN_Ip*q0*(fac2*z*g11+g12)
        gxz_(iz,eo) = dxdr*g12
        gyz_(iz,eo) = C_y*q0*SIGN_Ip*(fac2*z*g12+g22)
        gzz_(iz,eo) = g22

        ! jacobian
        jacobian_(iz,eo)=C_xy*abs(q0)*major_R*(1+trpeps*cost)**2

        ! Bfield
        Bfield_(iz,eo)   = B
        dBdx_(iz,eo)     = dBdr_c/dxdr
        dBdy_(iz,eo)     = 0._xp
        dBdz_(iz,eo)     = dBdchi

        !derivatives with respect to cylindrical coordinates
        dxdR_(iz,eo)  = cost
        dxdZ_(iz,eo)  = sint
        R_hat_(iz,eo) = major_R*(1.+trpeps*cost)
        Z_hat_(iz,eo) = major_Z+major_R*trpeps*sint

        !!!!!!! TEST salpha
    !         ! metric
    !   gxx_(iz,eo) = 1._xp
    !   gxy_(iz,eo) = shat*z - amhd*SIN(z)
    !   gxz_(iz,eo) = 0._xp
    !   gyy_(iz,eo) = 1._xp + (shat*z - amhd*SIN(z))**2
    !   gyz_(iz,eo) = 1._xp/trpeps
    !   gzz_(iz,eo) = 1._xp/trpeps**2
    !   dxdR_(iz,eo)= COS(z)
    !   dxdZ_(iz,eo)= SIN(z)

    ! ! ! Poloidal plane coordinates
    !   R_hat_(iz,eo) = 1._xp + trpeps*COS(z)
    !   Z_hat_(iz,eo) = trpeps*SIN(z)

    ! ! Relative strengh of modulus of B
    !   Bfield_(iz,eo) = 1._xp/(1._xp + trpeps*COS(z))

    ! ! Jacobian
    !   Jacobian_(iz,eo) = q0/Bfield_(iz,eo)

    ! ! Derivative of the magnetic field strenght
    !   dBdx_(iz,eo) = -COS(z)*Bfield_(iz,eo)**2 ! LB = 1
    !   dBdy_(iz,eo) =  0._xp
    !   dBdz_(iz,eo) =  trpeps*SIN(z)*Bfield_(iz,eo)**2

    ! ! Curvature factor
    !   C_xy = 1._xp
      !!!!! END TEST
      end do
    END DO

  END SUBROUTINE get_circ

  !-----------------------------------------------------------
END MODULE circular
