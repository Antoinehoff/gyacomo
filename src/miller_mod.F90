!! This source has been adapted from the GENE code https://genecode.org/ !!
!>Implementation of the local equilibrium model of [R.L. Miller et al., PoP 5, 973 (1998)
!>and [J. Candy, PPCF 51, 105009 (2009)]
MODULE miller
  USE prec_const
  USE basic
  USE parallel, ONLY: my_id, num_procs_z, nbr_U, nbr_D, comm0
  ! use coordinates,only: gcoor, get_dzprimedz
  USE grid, ONLY: local_Nky, local_Nkx, local_nz, ngz, nzgrid, kyarray, kxarray, zarray, total_nz, local_nz_offset, iodd, ieven, deltaz
  ! use discretization
  USE lagrange_interpolation

  implicit none
  public:: get_miller, set_miller_parameters
  public:: rho, kappa, delta, s_kappa, s_delta, drR, drZ, zeta, s_zeta
  public:: thetaShift
  public:: thetak, thetad

  private

  real(xp) :: rho, kappa, delta, s_kappa, s_delta, drR, drZ, zeta, s_zeta
  real(xp) :: thetaShift
  real(xp) :: thetak, thetad
  real(xp) :: asurf, theta1, theta2, theta3, delta2, delta3, Raxis, Zaxis
CONTAINS

  !>Set defaults for miller parameters
  subroutine set_miller_parameters(kappa_,s_kappa_,delta_,s_delta_,zeta_,s_zeta_,theta1_,theta2_)
    real(xp), INTENT(IN) :: kappa_,s_kappa_,delta_,s_delta_,zeta_,s_zeta_,theta1_,theta2_
    rho     = -1._xp
    kappa   = kappa_
    s_kappa = s_kappa_
    delta   = delta_
    s_delta = s_delta_
    zeta    = zeta_
    s_zeta  = s_zeta_
    drR     = 0._xp
    drZ     = 0._xp
    thetak  = 0._xp
    thetad  = 0._xp
    asurf   = 0.54_xp
    theta1  = theta1_
    theta2  = theta2_
    theta3  = 0._xp
    delta2  = 1._xp
    delta3  = 1._xp
    Raxis   = 1._xp
    Zaxis   = 0._xp

  end subroutine set_miller_parameters

  !>Get Miller metric, magnetic field, jacobian etc.
  subroutine get_miller(geom,trpeps,major_R,major_Z,q0,shat,Npol,amhd,edge_opt,&
       C_y,C_xy,Cyq0_x0,xpdx_pm_geom,gxx_,gxy_,gxz_,gyy_,gyz_,gzz_,dBdx_,dBdy_,dBdz_,&
       Bfield_,jacobian_,R_hat_,Z_hat_,dxdR_,dxdZ_)
    !!!!!!!!!!!!!!!! GYACOMO INTERFACE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER(len=16), INTENT(IN) :: geom
    real(xp), INTENT(INOUT) :: trpeps          ! eps in gyacomo (inverse aspect ratio)
    real(xp), INTENT(INOUT) :: major_R         ! major radius
    real(xp), INTENT(INOUT) :: major_Z         ! major Z
    real(xp), INTENT(INOUT) :: q0              ! safetyfactor
    real(xp), INTENT(INOUT) :: shat            ! safetyfactor
    INTEGER,  INTENT(IN)    :: Npol            ! number of poloidal turns
    real(xp), INTENT(INOUT) :: amhd            ! alpha mhd
    real(xp), INTENT(INOUT) :: edge_opt        ! optimization of point placement
    real(xp), INTENT(INOUT) :: xpdx_pm_geom    ! amplitude mag. eq. pressure grad.
    real(xp), INTENT(INOUT) :: C_y, C_xy, Cyq0_x0
    real(xp), dimension(local_nz+ngz,nzgrid), INTENT(OUT) :: &
                                              gxx_,gxy_,gxz_,gyy_,gyz_,gzz_,&
                                              dBdx_,dBdy_,dBdz_,Bfield_,jacobian_,&
                                              R_hat_,Z_hat_,dxdR_,dxdZ_
    INTEGER :: iz,eo
    ! No parameter in gyacomo yet
    real(xp) :: sign_Ip_CW=1._xp ! current sign (only normal current)
    real(xp) :: sign_Bt_CW=1._xp ! current sign (only normal current)
    !!!!!!!!!!!!!! END GYACOMO INTERFACE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer:: np, np_s, Npol_ext, Npol_s

    real(xp), dimension(500*(Npol+2)):: R,Z,R_rho,Z_rho,R_theta,Z_theta,R_theta_theta,Z_theta_theta,dlp,Rc,cosu,sinu,Bphi
    real(xp), dimension(500*(Npol+2)):: drRcirc, drRelong, drRelongTilt, drRtri, drRtriTilt, drZcirc, drZelong, drZelongTilt
    real(xp), dimension(500*(Npol+2)):: drZtri, drZtriTilt, dtdtRcirc, dtdtRelong, dtdtRelongTilt, dtdtRtri, dtdtRtriTilt
    real(xp), dimension(500*(Npol+2)):: dtdtZcirc, dtdtZelong, dtdtZelongTilt, dtdtZtri, dtdtZtriTilt, dtRcirc, dtRelong
    real(xp), dimension(500*(Npol+2)):: dtRelongTilt, dtRtri, dtRtriTilt, dtZcirc, dtZelong, dtZelongTilt, dtZtri, dtZtriTilt
    real(xp), dimension(500*(Npol+2)):: Rcirc, Relong, RelongTilt, Rtri, RtriTilt, Zcirc, Zelong, ZelongTilt, Ztri, ZtriTilt
    real(xp), dimension(500*(Npol+2)):: drrShape, drrAng, drxAng, dryAng, dtdtrShape, dtdtrAng, dtdtxAng
    real(xp), dimension(500*(Npol+2)):: dtdtyAng, dtrShape, dtrAng, dtxAng, dtyAng, rShape, rAng, xAng, yAng
    real(xp), dimension(500*(Npol+2)):: theta, thAdj, J_r, B, Bp, D0, D1, D2, D3, nu, chi, psi1, nu1
    real(xp), dimension(500*(Npol+2)):: tmp_reverse, theta_reverse, tmp_arr

    real(xp), dimension(500*(Npol+1)):: theta_s, thAdj_s, chi_s, theta_s_reverse
    real(xp), dimension(500*(Npol+1)):: R_s, Z_s, R_theta_s, Z_theta_s, Rc_s, cosu_s, sinu_s, Bphi_s, B_s, Bp_s, dlp_s
    real(xp), dimension(500*(Npol+1)):: dtRcirc_s, dtRelong_s, dtRelongTilt_s, dtRtri_s, dtRtriTilt_s, dtZcirc_s
    real(xp), dimension(500*(Npol+1)):: dtZelong_s, dtZelongTilt_s, dtZtri_s, dtZtriTilt_s, Rcirc_s, Relong_s, RelongTilt_s
    real(xp), dimension(500*(Npol+1)):: Rtri_s, RtriTilt_s, Zcirc_s, Zelong_s, ZelongTilt_s, Ztri_s, ZtriTilt_s!, dtrShape_s
    ! real(xp), dimension(500*(Npol+1)):: dtrAng_s, dtxAng_s, dtyAng_s, rShape_s, rAng_s, xAng_s, yAng_s
    real(xp), dimension(500*(Npol+1)):: psi1_s, nu1_s, dchidx_s, dB_drho_s, dB_dl_s, dnu_drho_s, dnu_dl_s, grad_nu_s
    real(xp), dimension(500*(Npol+1)):: gxx, gxy, gxz, gyy, gyz, gzz, dtheta_dchi_s, dBp_dchi_s, jacobian, dBdx, dBdz
    real(xp), dimension(500*(Npol+1)):: g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, tmp_arr_s, dxdR_s, dxdZ_s, K_x, K_y !tmp_arr2

    real(xp), dimension(1:total_nz+ngz):: gxx_out,gxy_out,gxz_out,gyy_out,gyz_out,gzz_out,Bfield_out,jacobian_out, dBdx_out, dBdz_out, chi_out
    real(xp), dimension(1:total_nz+ngz):: R_out, Z_out, dxdR_out, dxdZ_out
    real(xp):: d_inv, drPsi, dxPsi, dq_dx, dq_dpsi, R0, Z0, B0, F, D0_full, D1_full, D2_full, D3_full
    ! real(xp) :: Lnorm, Psi0 ! currently module-wide defined anyway
    real(xp):: pprime, ffprime, D0_mid, D1_mid, D2_mid, D3_mid, dx_drho, pi, mu_0, dzprimedz
    real(xp):: rho_a, psiN, drpsiN, CN2, CN3, Rcenter, Zcenter, drRcenter, drZcenter
    logical:: bMaxShift
    integer:: i, k, iBmax

    Npol_ext = Npol+2
    Npol_s   = Npol+1
    np   = 500*Npol_ext
    np_s = 500*Npol_s

    rho = trpeps*major_R
    if (rho.le.0._xp) ERROR STOP '>> ERROR << flux surface radius not defined'
    trpeps = rho/major_R

    q0 = sign_Ip_CW * sign_Bt_CW * abs(q0)

    R0=major_R
    B0=1._xp*sign_Bt_CW
    F=R0*B0
    Z0=major_Z
    pi = acos(-1._xp)
    mu_0=4._xp*pi

    theta=linspace(-pi*Npol_ext,pi*Npol_ext-2._xp*pi*Npol_ext/np,np)
    d_inv=asin(delta)

    thetaShift = 0._xp
    iBmax = 1
    bMaxShift = .true. ! limits the initialization routine to run no more than twice
    do while (bMaxShift)
      !flux surface parametrization
      thAdj = theta + thetaShift
      SELECT CASE (geom)
      CASE('Miller','miller')
        if (zeta/=0._xp .or. s_zeta/=0._xp) then
          R = R0 + rho*Cos(thAdj + d_inv*Sin(thAdj))
          Z = Z0 + kappa*rho*Sin(thAdj + zeta*Sin(2*thAdj))

          R_rho = drR + Cos(thAdj + d_inv*Sin(thAdj)) - s_delta*Sin(thAdj)*Sin(thAdj + d_inv*Sin(thAdj))
          Z_rho = drZ + kappa*s_zeta*Cos(thAdj + zeta*Sin(2*thAdj))*Sin(2*thAdj) &
                + kappa*Sin(thAdj + zeta*Sin(2*thAdj)) + kappa*s_kappa*Sin(thAdj + zeta*Sin(2*thAdj))

          R_theta = -(rho*(1 + d_inv*Cos(thAdj))*Sin(thAdj + d_inv*Sin(thAdj)))
          Z_theta = kappa*rho*(1 + 2*zeta*Cos(2*thAdj))*Cos(thAdj + zeta*Sin(2*thAdj))

          R_theta_theta = -(rho*(1 + d_inv*Cos(thAdj))**2*Cos(thAdj + d_inv*Sin(thAdj))) &
                + d_inv*rho*Sin(thAdj)*Sin(thAdj + d_inv*Sin(thAdj))
          Z_theta_theta = -4*kappa*rho*zeta*Cos(thAdj + zeta*Sin(2*thAdj))*Sin(2*thAdj) &
                - kappa*rho*(1 + 2*zeta*Cos(2*thAdj))**2*Sin(thAdj + zeta*Sin(2*thAdj))
        else
          Rcirc = rho*Cos(thAdj - thetad + thetak)
          Zcirc = rho*Sin(thAdj - thetad + thetak)
          Relong = Rcirc
          Zelong = Zcirc + (-1 + kappa)*rho*Sin(thAdj - thetad + thetak)
          RelongTilt = Relong*Cos(thetad - thetak) - Zelong*Sin(thetad - thetak)
          ZelongTilt = Zelong*Cos(thetad - thetak) + Relong*Sin(thetad - thetak)
          Rtri = RelongTilt - rho*Cos(thAdj) + rho*Cos(thAdj + delta*Sin(thAdj))
          Ztri = ZelongTilt
          RtriTilt = Rtri*Cos(thetad) + Ztri*Sin(thetad)
          ZtriTilt = Ztri*Cos(thetad) - Rtri*Sin(thetad)
          R = R0 + RtriTilt
          Z = Z0 + ZtriTilt

          drRcirc = Cos(thAdj - thetad + thetak)
          drZcirc = Sin(thAdj - thetad + thetak)
          drRelong = drRcirc
          drZelong = drZcirc - (1 - kappa - kappa*s_kappa)*Sin(thAdj - thetad + thetak)
          drRelongTilt = drRelong*Cos(thetad - thetak) - drZelong*Sin(thetad - thetak)
          drZelongTilt = drZelong*Cos(thetad - thetak) + drRelong*Sin(thetad - thetak)
          drRtri = drRelongTilt - Cos(thAdj) + Cos(thAdj + delta*Sin(thAdj)) &
                - s_delta*Sin(thAdj)*Sin(thAdj + delta*Sin(thAdj))
          drZtri = drZelongTilt
          drRtriTilt = drRtri*Cos(thetad) + drZtri*Sin(thetad)
          drZtriTilt = drZtri*Cos(thetad) - drRtri*Sin(thetad)
          R_rho = drR + drRtriTilt
          Z_rho = drZ + drZtriTilt

          dtRcirc = -(rho*Sin(thAdj - thetad + thetak))
          dtZcirc = rho*Cos(thAdj - thetad + thetak)
          dtRelong = dtRcirc
          dtZelong = dtZcirc + (-1 + kappa)*rho*Cos(thAdj - thetad + thetak)
          dtRelongTilt = dtRelong*Cos(thetad - thetak) - dtZelong*Sin(thetad - thetak)
          dtZelongTilt = dtZelong*Cos(thetad - thetak) + dtRelong*Sin(thetad - thetak)
          dtRtri = dtRelongTilt + rho*Sin(thAdj) - rho*(1 + delta*Cos(thAdj))*Sin(thAdj + delta*Sin(thAdj))
          dtZtri = dtZelongTilt
          dtRtriTilt = dtRtri*Cos(thetad) + dtZtri*Sin(thetad)
          dtZtriTilt = dtZtri*Cos(thetad) - dtRtri*Sin(thetad)
          R_theta = dtRtriTilt
          Z_theta = dtZtriTilt

          dtdtRcirc = -(rho*Cos(thAdj - thetad + thetak))
          dtdtZcirc = -(rho*Sin(thAdj - thetad + thetak))
          dtdtRelong = dtdtRcirc
          dtdtZelong = dtdtZcirc - (-1 + kappa)*rho*Sin(thAdj - thetad + thetak)
          dtdtRelongTilt = dtdtRelong*Cos(thetad - thetak) - dtdtZelong*Sin(thetad - thetak)
          dtdtZelongTilt = dtdtZelong*Cos(thetad - thetak) + dtdtRelong*Sin(thetad - thetak)
          dtdtRtri = dtdtRelongTilt + rho*Cos(thAdj) - rho*(1 + delta*Cos(thAdj))**2*Cos(thAdj + delta*Sin(thAdj)) &
                + delta*rho*Sin(thAdj)*Sin(thAdj + delta*Sin(thAdj))
          dtdtZtri = dtdtZelongTilt
          dtdtRtriTilt = dtdtRtri*Cos(thetad) + dtdtZtri*Sin(thetad)
          dtdtZtriTilt = dtdtZtri*Cos(thetad) - dtdtRtri*Sin(thetad)
          R_theta_theta = dtdtRtriTilt
          Z_theta_theta = dtdtZtriTilt
        endif
      CASE('Miller_global','miller_global')
        rho_a = rho/aSurf

        CN2 = (-1._xp + Delta2**2)/(1._xp + Delta2**2)
        CN3 = (-1._xp + Delta3**2)/(1._xp + Delta3**3)

        psiN = rho_a**2 + CN2*rho_a**2 + CN3*rho_a**3
        drpsiN = 2._xp*rho_a + 2._xp*CN2*rho_a + 3._xp*CN3*rho_a**2

        xAng = 2*(1._xp + CN2*Cos(2._xp*(thAdj + theta2)))**3 - 27._xp*CN3**2*Cos(3._xp*(thAdj + theta3))**2*psiN
        drxAng = -27._xp*CN3**2*Cos(3._xp*(thAdj + theta3))**2*drpsiN
        dtxAng = -12._xp*CN2*(1._xp + CN2*Cos(2._xp*(thAdj + theta2)))**2*Sin(2._xp*(thAdj + theta2)) &
             + 162._xp*CN3**2*Cos(3._xp*(thAdj + theta3))*psiN*Sin(3._xp*(thAdj + theta3))
        dtdtxAng = 486._xp*CN3**2*Cos(6._xp*(thAdj + theta3))*psiN + 24._xp*CN2*(1._xp + CN2*Cos(2._xp*(thAdj + theta2))) &
             *(-(Cos(2._xp*(thAdj + theta2))*(1._xp + CN2*Cos(2._xp*(thAdj + theta2)))) + 2._xp*CN2*Sin(2._xp*(thAdj + theta2))**2)

        yAng = 3._xp*Sqrt(3._xp)*CN3*Cos(3._xp*(thAdj + theta3))*Sqrt(psiN)*Sqrt(2._xp*(1._xp + CN2*Cos(2._xp*(thAdj + theta2)))**3 + xAng)
        dryAng = (yAng*(drpsiN/psiN + drxAng/(2._xp*(1._xp + CN2*Cos(2._xp*(thAdj + theta2)))**3 + xAng)))/2._xp
        dtyAng = yAng*(-3._xp*Tan(3._xp*(thAdj + theta3)) &
             + (-12._xp*CN2*(1._xp + CN2*Cos(2._xp*(thAdj + theta2)))**2*Sin(2._xp*(thAdj + theta2)) + dtxAng) &
             /(2._xp*(2*(1._xp + CN2*Cos(2*(thAdj + theta2)))**3 + xAng)))
        dtdtyAng = dtyAng**2/yAng + yAng*(-9._xp/Cos(3._xp*(thAdj + theta3))**2 &
             - (-12._xp*CN2*(1._xp + CN2*Cos(2._xp*(thAdj + theta2)))**2*Sin(2._xp*(thAdj + theta2)) + dtxAng)**2 &
             /(2._xp*(2._xp*(1._xp + CN2*Cos(2._xp*(thAdj + theta2)))**3 + xAng)**2) &
             + (-24._xp*CN2*Cos(2._xp*(thAdj + theta2))*(1._xp + CN2*Cos(2._xp*(thAdj + theta2)))**2 &
             + 48._xp*CN2**2*(1._xp + CN2*Cos(2._xp*(thAdj + theta2)))*Sin(2._xp*(thAdj + theta2))**2 + dtdtxAng) &
             /(2._xp*(2._xp*(1._xp + CN2*Cos(2._xp*(thAdj + theta2)))**3 + xAng)))

        rAng = atan2(yAng,xAng)/3.
        drrAng = (-(yAng*drxAng) + xAng*dryAng)/(3._xp*(xAng**2 + yAng**2))
        dtrAng = (-(yAng*dtxAng) + xAng*dtyAng)/(3._xp*(xAng**2 + yAng**2))
        dtdtrAng = (-6._xp*dtrAng**2*yAng)/xAng &
             + ((2._xp*yAng*dtxAng**2)/xAng - 2._xp*dtxAng*dtyAng - yAng*dtdtxAng + xAng*dtdtyAng)/(3._xp*(xAng**2 + yAng**2))

        rShape = (aSurf*(1._xp + CN2*Cos(2._xp*(thAdj + theta2)))/Cos(3._xp*(thAdj + theta3))*&
             &(-1._xp + Cos(rAng) + Sqrt(3._xp)*Sin(rAng)))/(3._xp*CN3)
        drrShape = (aSurf*(1._xp + CN2*Cos(2._xp*(thAdj + theta2)))/Cos(3._xp*(thAdj + theta3))*&
             &(Sqrt(3._xp)*Cos(rAng) - Sin(rAng))*drrAng)/(3._xp*CN3)
        dtrShape = rShape*((-2._xp*CN2*Sin(2._xp*(thAdj + theta2)))/(1._xp + CN2*Cos(2._xp*(thAdj + theta2)))&
             &+ 3._xp*Tan(3._xp*(thAdj + theta3)) &
             &+ ((Sqrt(3._xp)*Cos(rAng) - Sin(rAng))*dtrAng)/(-1._xp + Cos(rAng) + Sqrt(3._xp)*Sin(rAng)))
        dtdtrShape = dtrShape**2/rShape + rShape*(9._xp - (4._xp*CN2*Cos(2._xp*(thAdj + theta2)))/(1._xp + CN2*Cos(2._xp*(thAdj + theta2))) &
             - (4._xp*CN2**2*Sin (2._xp*(thAdj + theta2))**2)/(1._xp + CN2*Cos(2._xp*(thAdj + theta2)))**2 + 9._xp*Tan (3._xp*(thAdj + theta3))**2 &
             + ((-4._xp + Cos(rAng) + Sqrt(3._xp)*Sin(rAng))*dtrAng**2)/(-1._xp + Cos(rAng) + Sqrt(3._xp)*Sin(rAng))**2 &
             + ((Sqrt(3._xp)*Cos(rAng) - Sin(rAng))*dtdtrAng)/(-1._xp + Cos(rAng) + Sqrt(3._xp)*Sin(rAng)))

        do i=1._xp,np
           if (abs(CN3*cos(3._xp*(thAdj(i)+theta3)))<1e-8) then
              rShape(i) = aSurf*Sqrt(psiN/(1._xp + CN2*Cos(2*(thAdj(i) + theta2))))
              drrShape(i) = (rShape(i)*drpsiN)/(2._xp*psiN)
              dtrShape(i) = (CN2*rShape(i)*Sin(2._xp*(thAdj(i) + theta2)))/(1._xp + CN2*Cos(2._xp*(thAdj(i) + theta2)))
              dtdtrShape(i) = dtrShape(i)**2/rShape(i) &
                   + (2*(CN2**2 + CN2*Cos(2._xp*(thAdj(i) + theta2)))*rShape(i))/(1._xp + CN2*Cos(2._xp*(thAdj(i) + theta2)))**2
           endif

           if (abs(1._xp + CN2*Cos(2._xp*(thAdj(i) + theta2)))<1e-8) then
              rShape(i) = aSurf*Sqrt(psiN)
              drrShape(i) = 1._xp
              dtrShape(i) = 0._xp
              dtdtrShape(i) = 0._xp
           endif
        enddo

        Rcenter = Raxis - (-R0 + Raxis)*rho_a**2
        Zcenter = Zaxis - rho_a**2*(-Z0 + Zaxis)
        drRcenter = -2._xp*(-R0 + Raxis)*rho_a
        drZcenter = -2._xp*rho_a*(-Z0 + Zaxis)

        R = Rcenter + Cos(thAdj)*rShape
        Z = rShape*Sin(thAdj) + Zcenter
        R_rho = (drRcenter + Cos(thAdj)*drrShape)/aSurf ! Adjust for rho deriv, not rho_a deriv
        Z_rho = (drZcenter + Sin(thAdj)*drrShape)/aSurf ! Adjust for rho deriv, not rho_a deriv
        R_theta = -(rShape*Sin(thAdj)) + Cos(thAdj)*dtrShape
        Z_theta = Cos(thAdj)*rShape + Sin(thAdj)*dtrShape
        R_theta_theta = -(Cos(thAdj)*rShape) - 2._xp*Sin(thAdj)*dtrShape + Cos(thAdj)*dtdtrShape
        Z_theta_theta = -(rShape*Sin(thAdj)) + 2._xp*Cos(thAdj)*dtrShape + Sin(thAdj)*dtdtrShape

      END SELECT

      !dl/dtheta
      dlp=(R_theta**2+Z_theta**2)**0.5_xp

      !curvature radius
      Rc=dlp**3*(R_theta*Z_theta_theta-Z_theta*R_theta_theta)**(-1)

      ! some useful quantities (see papers for definition of u)
      cosu=Z_theta/dlp
      sinu=-R_theta/dlp

      !Jacobian J_r = (dPsi/dr) J_psi = (dPsi/dr) / [(nabla fz x nabla psi)* nabla theta]
      !             = R * (dR/drho dZ/dtheta - dR/dtheta dZ/drho) = R dlp / |nabla r|
      J_r=R*(R_rho*Z_theta-R_theta*Z_rho)

      !From definition of q = 1/(2 pi) int (B nabla fz) / (B nabla theta) dtheta:
      !dPsi/dr = sign_Bt sign_Ip / (2 pi q) int F / R^2 J_r dtheta
      !        = F / (2 pi |q|) int J_r/R^2 dtheta
      tmp_arr=J_r/R**2
      drPsi=sign_Ip_CW*F/(2._xp*pi*Npol_ext*q0)*sum(tmp_arr)*2._xp*pi*Npol_ext/np !dlp_int(tmp_arr,1.0)

      !Poloidal field (Bp = Bvec * nabla l)
      Bp=sign_Ip_CW * drPsi / J_r * dlp

      !toroidal field
      Bphi=F/R

      !total modulus of Bfield
      B=sqrt(Bphi**2+Bp**2)

      bMaxShift = .false.
      if (thetaShift==0._xp) then
        do i = 2,np-1
          if (B(iBmax)<B(i)) then
              iBmax = i
          end if
        enddo
        if (iBmax/=1) then
          bMaxShift = .true.
          thetaShift = theta(iBmax)-theta(1)
        end if
      end if
    enddo

    !definition of radial coordinate! dx_drho=1 --> x = r
    dx_drho=1. !drPsi/Psi0*Lnorm*q0
    CALL speak('Using radial coordinate with dx/dr = '//str(dx_drho),2)
    dxPsi  = drPsi/dx_drho
    C_y    = dxPsi*sign_Ip_CW
    C_xy   = abs(B0*dxPsi/C_y)

    CALL speak("Setting C_xy = "//str(C_xy)//' C_y = '//str(C_y)//" C_x' = "//str(1./dxPsi),2)
    CALL speak("B_unit/Bref conversion factor = "//str(q0/rho*drPsi),2)
    CALL speak("dPsi/dr = "//str(drPsi),2)
    IF (thetaShift.ne.0._xp) CALL speak("thetaShift = "//str(thetaShift),2)


    !--------shear is expected to be defined as rho/q*dq/drho--------!
    dq_dx   = shat*q0/rho/dx_drho
    Cyq0_x0 = C_y*q0/rho
    dq_dpsi = dq_dx/dxPsi
    pprime  = -amhd/q0**2/R0/(2*mu_0)*B0**2/drPsi

    !neg. xpdx normalized to magnetic pressure for pressure term
    xpdx_pm_geom=amhd/q0**2/R0/dx_drho

    !first coefficient of psi in varrho expansion
    psi1 = R*Bp*sign_Ip_CW

    !integrals for ffprime evaluation
    do i=1,np
       tmp_arr=(2._xp/Rc-2._xp*cosu/R)/(R*psi1**2)
       D0(i)=-F*dlp_int_ind(tmp_arr,dlp,i)
       tmp_arr=B**2*R/psi1**3
       D1(i)=-dlp_int_ind(tmp_arr,dlp,i)/F
       tmp_arr=mu_0*R/psi1**3
       D2(i)=-dlp_int_ind(tmp_arr,dlp,i)*F
       tmp_arr=1./(R*psi1)
       D3(i)=-dlp_int_ind(tmp_arr,dlp,i)*F
    enddo
    tmp_arr=(2._xp/Rc-2._xp*cosu/R)/(R*psi1**2)
    D0_full=-F*dlp_int(tmp_arr,dlp)
    tmp_arr=B**2*R/psi1**3
    D1_full=-dlp_int(tmp_arr,dlp)/F
    tmp_arr=mu_0*R/psi1**3
    D2_full=-dlp_int(tmp_arr,dlp)*F
    tmp_arr=1./(R*psi1)
    D3_full=-dlp_int(tmp_arr,dlp)*F
    D0_mid=D0(np/2+1)
    D1_mid=D1(np/2+1)
    D2_mid=D2(np/2+1)
    D3_mid=D3(np/2+1)

    ffprime=-(sign_Ip_CW*dq_dpsi*2._xp*pi*Npol_ext+D0_full+D2_full*pprime)/D1_full

    CALL speak("ffprime = "//str(ffprime),2)

    D0=D0-D0_mid
    D1=D1-D1_mid
    D2=D2-D2_mid
    nu=D3-D3_mid

    nu1=psi1*(D0+D1*ffprime+D2*pprime)

    !straight field line angle defined on equidistant theta grid
    !alpha = fz + nu = - (q chi - fz) => chi = -nu / q
    chi=-nu/q0

    !correct small scaling error (<0.5%, due to finite integration resolution)
    chi=chi*(maxval(theta)-minval(theta))/(maxval(chi)-minval(chi))

    !new grid equidistant in straight field line angle
    chi_s = linspace(-pi*Npol_s,pi*Npol_s-2*pi*Npol_s/np_s,np_s)

    if (sign_Ip_CW.lt.0._xp) then !make chi increasing function to not confuse lag3interp
       tmp_reverse = chi(np:1:-1)
       theta_reverse = theta(np:1:-1)
       call lag3interp(theta_reverse,tmp_reverse,np,theta_s,chi_s,np_s)
       theta_s_reverse = theta_s(np_s:1:-1)
    else
       !lag3interp(y_in,x_in,n_in,y_out,x_out,n_out)
       call lag3interp(theta,chi,np,theta_s,chi_s,np_s)
    endif
    dtheta_dchi_s=deriv_fd(theta_s,chi_s,np_s)

    !arrays equidistant in straight field line angle
    thAdj_s = theta_s + thetaShift

    if (zeta/=0._xp .or. s_zeta/=0._xp) then
      R_s = R0 + rho*Cos(thAdj_s + d_inv*Sin(thAdj_s))
      Z_s = Z0 + kappa*rho*Sin(thAdj_s + zeta*Sin(2*thAdj_s))

      R_theta_s = -(dtheta_dchi_s*rho*(1 + d_inv*Cos(thAdj_s))*Sin(thAdj_s + d_inv*Sin(thAdj_s)))
      Z_theta_s = dtheta_dchi_s*kappa*rho*(1 + 2*zeta*Cos(2*thAdj_s))*Cos(thAdj_s + zeta*Sin(2*thAdj_s))
    else
      Rcirc_s = rho*Cos(thAdj_s - thetad + thetak)
      Zcirc_s = rho*Sin(thAdj_s - thetad + thetak)
      Relong_s = Rcirc_s
      Zelong_s = Zcirc_s + (-1 + kappa)*rho*Sin(thAdj_s - thetad + thetak)
      RelongTilt_s = Relong_s*Cos(thetad - thetak) - Zelong_s*Sin(thetad - thetak)
      ZelongTilt_s = Zelong_s*Cos(thetad - thetak) + Relong_s*Sin(thetad - thetak)
      Rtri_s = RelongTilt_s - rho*Cos(thAdj_s) + rho*Cos(thAdj_s + delta*Sin(thAdj_s))
      Ztri_s = ZelongTilt_s
      RtriTilt_s = Rtri_s*Cos(thetad) + Ztri_s*Sin(thetad)
      ZtriTilt_s = Ztri_s*Cos(thetad) - Rtri_s*Sin(thetad)
      R_s = R0 + RtriTilt_s
      Z_s = Z0 + ZtriTilt_s

      dtRcirc_s = -(rho*Sin(thAdj_s - thetad + thetak))
      dtZcirc_s =   rho*Cos(thAdj_s - thetad + thetak)
      dtRelong_s = dtRcirc_s
      dtZelong_s = dtZcirc_s + (-1 + kappa)*rho*Cos(thAdj_s - thetad + thetak)
      dtRelongTilt_s = dtRelong_s*Cos(thetad - thetak) - dtZelong_s*Sin(thetad - thetak)
      dtZelongTilt_s = dtZelong_s*Cos(thetad - thetak) + dtRelong_s*Sin(thetad - thetak)
      dtRtri_s = dtRelongTilt_s + rho*Sin(thAdj_s) &
           - rho*(1 + delta*Cos(thAdj_s))*Sin(thAdj_s + delta*Sin(thAdj_s))
      dtZtri_s = dtZelongTilt_s
      dtRtriTilt_s = dtRtri_s*Cos(thetad) + dtZtri_s*Sin(thetad)
      dtZtriTilt_s = dtZtri_s*Cos(thetad) - dtRtri_s*Sin(thetad)
      R_theta_s = dtheta_dchi_s*dtRtriTilt_s
      Z_theta_s = dtheta_dchi_s*dtZtriTilt_s
    endif
    if (sign_Ip_CW.lt.0._xp) then
       call lag3interp(nu1,theta,np,tmp_arr_s,theta_s_reverse,np_s)
       nu1_s = tmp_arr_s(np_s:1:-1)
       call lag3interp(Bp,theta,np,tmp_arr_s,theta_s_reverse,np_s)
       Bp_s = tmp_arr_s(np_s:1:-1)
       call lag3interp(dlp,theta,np,tmp_arr_s,theta_s_reverse,np_s)
       dlp_s = tmp_arr_s(np_s:1:-1)
       call lag3interp(Rc,theta,np,tmp_arr_s,theta_s_reverse,np_s)
       Rc_s = tmp_arr_s(np_s:1:-1)
    else
       call lag3interp(nu1,theta,np,nu1_s,theta_s,np_s)
       call lag3interp(Bp,theta,np,Bp_s,theta_s,np_s)
       call lag3interp(dlp,theta,np,dlp_s,theta_s,np_s)
       call lag3interp(Rc,theta,np,Rc_s,theta_s,np_s)
    endif

    psi1_s = R_s*Bp_s*sign_Ip_CW

    dBp_dchi_s=deriv_fd(Bp_s,chi_s,np_s)

    Bphi_s=F/R_s
    B_s=sqrt(Bphi_s**2+Bp_s**2)
    cosu_s=Z_theta_s/dlp_s/dtheta_dchi_s
    sinu_s=-R_theta_s/dlp_s/dtheta_dchi_s

    !radial derivative of straight field line angle
    dchidx_s=-(nu1_s/psi1_s*dxPsi+chi_s*dq_dx)/q0

    !Bfield derivatives in Mercier-Luc coordinates (varrho,l,fz)
    dB_drho_s=-1./B_s*(F**2/R_s**3*cosu_s+Bp_s**2/Rc_s+mu_0*psi1_s*pprime)
    dB_dl_s=1./B_s*(Bp_s*dBp_dchi_s/dtheta_dchi_s/dlp_s+F**2/R_s**3*sinu_s)

    dnu_drho_s=nu1_s
    dnu_dl_s=-F/(R_s*psi1_s)
    grad_nu_s=sqrt(dnu_drho_s**2+dnu_dl_s**2)

    !contravariant metric coefficients (varrho,l,phi)->(x,y,z)
    gxx=(psi1_s/dxPsi)**2
    gxy=-psi1_s/dxPsi*C_y*sign_Ip_CW*nu1_s
    gxz=-psi1_s/dxPsi*(nu1_s+psi1_s*dq_dpsi*chi_s)/q0
    gyy=C_y**2*(grad_nu_s**2+1/R_s**2)
    gyz=sign_Ip_CW*C_y/q0*(grad_nu_s**2+dq_dpsi*nu1_s*psi1_s*chi_s)
    gzz=1./q0**2*(grad_nu_s**2+2.*dq_dpsi*nu1_s*psi1_s*chi_s+(dq_dpsi*psi1_s*chi_s)**2)

    jacobian=1./sqrt(gxx*gyy*gzz + 2.*gxy*gyz*gxz - gxz**2*gyy - gyz**2*gxx - gzz*gxy**2)

    !covariant metric coefficients
    g_xx=jacobian**2*(gyy*gzz-gyz**2)
    g_xy=jacobian**2*(gxz*gyz-gxy*gzz)
    g_xz=jacobian**2*(gxy*gyz-gxz*gyy)
    g_yy=jacobian**2*(gxx*gzz-gxz**2)
    g_yz=jacobian**2*(gxz*gxy-gxx*gyz)
    g_zz=jacobian**2*(gxx*gyy-gxy**2)

    !Bfield derivatives
    !dBdx = e_x * nabla B = J (nabla y x nabla z) * nabla B
    dBdx=jacobian*C_y/(q0*R_s)*(F/(R_s*psi1_s)*dB_drho_s+(nu1_s+dq_dpsi*chi_s*psi1_s)*dB_dl_s)
    dBdz=1./B_s*(Bp_s*dBp_dchi_s-F**2/R_s**3*R_theta_s)

    !curvature terms (these are just local and will be recalculated in geometry module)
    K_x = (0._xp-g_yz/g_zz*dBdz)
    K_y = (dBdx-g_xz/g_zz*dBdz)

    !(R,Z) derivatives for visualization
    dxdR_s = dx_drho/drPsi*psi1_s*cosu_s
    dxdZ_s = dx_drho/drPsi*psi1_s*sinu_s

    ! GHOSTS ADAPTED VERSION
    if (edge_opt==0._xp) then
      !gyacomo z-grid wo ghosts
      ! chi_out=linspace(-pi*Npol,pi*Npol-2*pi*Npol/total_nz,total_nz)
      !gyacomo z-grid with ghosts
      chi_out=linspace(-pi*Npol-4._xp*pi*Npol/total_nz,pi*Npol+2._xp*pi*Npol/total_nz,total_nz+ngz)
    else
      ERROR STOP '>> ERROR << ghosts not implemented for edge_opt yet'
      !new parallel coordinate chi_out==zprime
      !see also tracer_aux.F90
      if (Npol>1) ERROR STOP '>> ERROR << Npol>1 has not been implemented for edge_opt=\=0._xp'
      do k=1,total_nz
          chi_out(k)=sinh((-pi+k*2._xp*pi/total_nz)*log(edge_opt*pi+sqrt(edge_opt**2*pi**2+1))/pi)/edge_opt
      enddo
      !transform metrics according to chain rule
      do k=1,np_s
        !>dz'/dz conversion for edge_opt as function of z
          if (edge_opt.gt.0) then
            dzprimedz = edge_opt*pi/log(edge_opt*pi+sqrt((edge_opt*pi)**2+1._xp))/&
                  sqrt((edge_opt*chi_s(k))**2+1)
          else
            dzprimedz = 1.0
          endif
          gxz(k)=gxz(k)*dzprimedz
          gyz(k)=gyz(k)*dzprimedz
          gzz(k)=gzz(k)*dzprimedz**2
          jacobian(k)=jacobian(k)/dzprimedz
          dBdz(k)=dBdz(k)/dzprimedz
      enddo
    endif !edge_opt

    !! Loop over the z-grids and interpolate the results on it
    do eo=ieven,iodd
      ! Shift the grid
      chi_out = chi_out + REAL(eo-1,xp)*0.5*deltaz
      ! interpolate with ghosts
      call lag3interp(gxx,chi_s,np_s,gxx_out,chi_out,total_nz+ngz)
      call lag3interp(gxy,chi_s,np_s,gxy_out,chi_out,total_nz+ngz)
      call lag3interp(gxz,chi_s,np_s,gxz_out,chi_out,total_nz+ngz)
      call lag3interp(gyy,chi_s,np_s,gyy_out,chi_out,total_nz+ngz)
      call lag3interp(gyz,chi_s,np_s,gyz_out,chi_out,total_nz+ngz)
      call lag3interp(gzz,chi_s,np_s,gzz_out,chi_out,total_nz+ngz)
      call lag3interp(B_s,chi_s,np_s,Bfield_out,chi_out,total_nz+ngz)
      call lag3interp(dBdx,chi_s,np_s,dBdx_out,chi_out,total_nz+ngz)
      call lag3interp(dBdz,chi_s,np_s,dBdz_out,chi_out,total_nz+ngz)
      call lag3interp(jacobian,chi_s,np_s,jacobian_out,chi_out,total_nz+ngz)
      call lag3interp(R_s,chi_s,np_s,R_out,chi_out,total_nz+ngz)
      call lag3interp(Z_s,chi_s,np_s,Z_out,chi_out,total_nz+ngz)
      call lag3interp(dxdR_s,chi_s,np_s,dxdR_out,chi_out,total_nz+ngz)
      call lag3interp(dxdZ_s,chi_s,np_s,dxdZ_out,chi_out,total_nz+ngz)
      ! Fill the local geom arrays with the results
      DO iz = 1,local_nz + ngz
        gxx_(iz,eo)      = gxx_out(iz+local_nz_offset)
        gxy_(iz,eo)      = gxy_out(iz+local_nz_offset)
        gxz_(iz,eo)      = gxz_out(iz+local_nz_offset)
        gyy_(iz,eo)      = gyy_out(iz+local_nz_offset)
        gyz_(iz,eo)      = gyz_out(iz+local_nz_offset)
        gzz_(iz,eo)      = gzz_out(iz+local_nz_offset)
        Bfield_(iz,eo)   = Bfield_out(iz+local_nz_offset)
        dBdx_(iz,eo)     = dBdx_out(iz+local_nz_offset)
        dBdy_(iz,eo)     = 0._xp
        dBdz_(iz,eo)     = dBdz_out(iz+local_nz_offset)
        jacobian_(iz,eo) = jacobian_out(iz+local_nz_offset)
        R_hat_(iz,eo)    = R_out(iz+local_nz_offset)
        Z_hat_(iz,eo)    = Z_out(iz+local_nz_offset)
        dxdR_(iz,eo)     = dxdR_out(iz+local_nz_offset)
        dxdZ_(iz,eo)     = dxdZ_out(iz+local_nz_offset)
      ENDDO
  ENDDO
  contains
    !> Generate an equidistant array from min to max with n points
    function linspace(min,max,n) result(out)
      real(xp), INTENT(IN):: min, max
      integer,  INTENT(IN):: n
      real(xp), dimension(n):: out

      do i=1,n
         out(i)=min+(i-1)*(max-min)/(n-1)
      enddo
    end function linspace

    !> Weighted average
    real(xp) function average(var,weight)
      real(xp), dimension(np), INTENT(IN):: var, weight
      average=sum(var*weight)/sum(weight)
    end function average

    !> full theta integral with weight function dlp
    real(xp) function dlp_int(var,dlp)
      real(xp), dimension(np), INTENT(IN):: var, dlp
      dlp_int=sum(var*dlp)*2._xp*pi*Npol_ext/np
    end function dlp_int

    !> theta integral with weight function dlp, up to index 'ind'
    real(xp) function dlp_int_ind(var,dlp,ind)
      real(xp), dimension(np), INTENT(IN):: var, dlp
      integer, INTENT(IN):: ind
      dlp_int_ind=0._xp
      if (ind.gt.1) then
         dlp_int_ind=dlp_int_ind+var(1)*dlp(1)*pi*Npol_ext/np
         dlp_int_ind=dlp_int_ind+(sum(var(2:ind-1)*dlp(2:ind-1)))*2*pi*Npol_ext/np
         dlp_int_ind=dlp_int_ind+var(ind)*dlp(ind)*pi*Npol_ext/np
      endif
    end function dlp_int_ind

    !> 1st derivative with 2nd order finite differences
    function deriv_fd(y,x,n) result(out)
      integer,                INTENT(IN) :: n
      real(xp), dimension(n), INTENT(IN):: x,y
      real(xp), dimension(n) :: out,dx
      !call lag3deriv(y,x,n,out,x,n)
      out=0._xp
      do i=2,n-1
         out(i)=out(i)-y(i-1)/2._xp
         out(i)=out(i)+y(i+1)/2._xp
      enddo
      out(1)=y(2)-y(1)
      out(n)=y(n)-y(n-1)
      dx=x(2)-x(1)
      out=out/dx
    end function deriv_fd

  end subroutine get_miller


END MODULE miller
