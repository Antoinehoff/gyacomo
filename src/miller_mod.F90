!! This source has been adapted from GENE https://genecode.org/ !!
!>Implementation of the local equilibrium model of [R.L. Miller et al., PoP 5, 973 (1998)
!>and [J. Candy, PPCF 51, 105009 (2009)]
MODULE miller
  USE prec_const
  USE basic
  ! use coordinates,only: gcoor, get_dzprimedz
  USE grid
  ! use discretization
  USE lagrange_interpolation
  ! use par_in, only: beta, sign_Ip_CW, sign_Bt_CW, Npol
  USE model

  implicit none
  public:: get_miller, set_miller_parameters
  public:: rho, kappa, delta, s_kappa, s_delta, drR, drZ, zeta, s_zeta
  public:: thetaShift
  public:: thetak, thetad

  private

  real(dp) :: rho, kappa, delta, s_kappa, s_delta, drR, drZ, zeta, s_zeta
  real(dp) :: thetaShift
  real(dp) :: thetak, thetad

CONTAINS

  !>Set defaults for miller parameters
  subroutine set_miller_parameters(kappa_,s_kappa_,delta_,s_delta_,zeta_,s_zeta_)
    real(dp), INTENT(IN) :: kappa_,s_kappa_,delta_,s_delta_,zeta_,s_zeta_
    rho     = -1.0
    kappa   = kappa_
    s_kappa = s_kappa_
    delta   = delta_
    s_delta = s_delta_
    zeta    = zeta_
    s_zeta  = s_zeta_
    drR     = 0.0
    drZ     = 0.0

    thetak = 0.0
    thetad = 0.0

  end subroutine set_miller_parameters

  !>Get Miller metric, magnetic field, jacobian etc.
  subroutine get_miller(trpeps,major_R,major_Z,q0,shat,amhd,edge_opt,&
       C_y,C_xy,dpdx_pm_geom,gxx_,gyy_,gzz_,gxy_,gxz_,gyz_,dBdx_,dBdy_,&
       Bfield_,jacobian_,dBdz_,R_hat_,Z_hat_,dxdR_,dxdZ_,Ckxky_,gradz_coeff_)
    !!!!!!!!!!!!!!!! GYACOMO INTERFACE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(dp), INTENT(INOUT) :: trpeps          ! eps in gyacomo (inverse aspect ratio)
    real(dp), INTENT(INOUT) :: major_R         ! major radius
    real(dp), INTENT(INOUT) :: major_Z         ! major Z
    real(dp), INTENT(INOUT) :: q0              ! safetyfactor
    real(dp), INTENT(INOUT) :: shat            ! safetyfactor
    real(dp), INTENT(INOUT) :: amhd            ! alpha mhd
    real(dp), INTENT(INOUT) :: edge_opt        ! alpha mhd
    real(dp), INTENT(INOUT) :: dpdx_pm_geom    ! amplitude mag. eq. pressure grad.
    real(dp), INTENT(INOUT) :: C_y, C_xy

    real(dp), dimension(izgs:izge,0:1), INTENT(INOUT) :: &
                                              gxx_,gyy_,gzz_,gxy_,gxz_,gyz_,&
                                              dBdx_,dBdy_,Bfield_,jacobian_,&
                                              dBdz_,R_hat_,Z_hat_,dxdR_,dxdZ_, &
                                              gradz_coeff_
    real(dp), dimension(ikys:ikye,ikxs:ikxe,izgs:izge,0:1), INTENT(INOUT) :: Ckxky_
    ! No parameter in gyacomo yet
    real(dp) :: sign_Ip_CW=1 ! current sign (only normal current)
    real(dp) :: sign_Bt_CW=1 ! current sign (only normal current)
    !!!!!!!!!!!!!! END GYACOMO INTERFACE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Auxiliary variables for curvature computation
    real(dp) :: G1,G2,G3,Cx,Cy,ky,kx

    integer:: np, np_s, Npol_ext, Npol_s

    real(dp), dimension(500*(Npol+2)):: R,Z,R_rho,Z_rho,R_theta,Z_theta,R_theta_theta,Z_theta_theta,dlp,Rc,cosu,sinu,Bphi
    real(dp), dimension(500*(Npol+2)):: drRcirc, drRelong, drRelongTilt, drRtri, drRtriTilt, drZcirc, drZelong, drZelongTilt
    real(dp), dimension(500*(Npol+2)):: drZtri, drZtriTilt, dtdtRcirc, dtdtRelong, dtdtRelongTilt, dtdtRtri, dtdtRtriTilt
    real(dp), dimension(500*(Npol+2)):: dtdtZcirc, dtdtZelong, dtdtZelongTilt, dtdtZtri, dtdtZtriTilt, dtRcirc, dtRelong
    real(dp), dimension(500*(Npol+2)):: dtRelongTilt, dtRtri, dtRtriTilt, dtZcirc, dtZelong, dtZelongTilt, dtZtri, dtZtriTilt
    real(dp), dimension(500*(Npol+2)):: Rcirc, Relong, RelongTilt, Rtri, RtriTilt, Zcirc, Zelong, ZelongTilt, Ztri, ZtriTilt
    ! real(dp), dimension(500*(Npol+2)):: drrShape, drrAng, drxAng, dryAng, dtdtrShape, dtdtrAng, dtdtxAng
    ! real(dp), dimension(500*(Npol+2)):: dtdtyAng, dtrShape, dtrAng, dtxAng, dtyAng, rShape, rAng, xAng, yAng
    real(dp), dimension(500*(Npol+2)):: theta, thAdj, J_r, B, Bp, D0, D1, D2, D3, nu, chi, psi1, nu1
    real(dp), dimension(500*(Npol+2)):: tmp_reverse, theta_reverse, tmp_arr

    real(dp), dimension(500*(Npol+1)):: theta_s, thAdj_s, chi_s, theta_s_reverse
    real(dp), dimension(500*(Npol+1)):: R_s, Z_s, R_theta_s, Z_theta_s, Rc_s, cosu_s, sinu_s, Bphi_s, B_s, Bp_s, dlp_s
    real(dp), dimension(500*(Npol+1)):: dtRcirc_s, dtRelong_s, dtRelongTilt_s, dtRtri_s, dtRtriTilt_s, dtZcirc_s
    real(dp), dimension(500*(Npol+1)):: dtZelong_s, dtZelongTilt_s, dtZtri_s, dtZtriTilt_s, Rcirc_s, Relong_s, RelongTilt_s
    real(dp), dimension(500*(Npol+1)):: Rtri_s, RtriTilt_s, Zcirc_s, Zelong_s, ZelongTilt_s, Ztri_s, ZtriTilt_s!, dtrShape_s
    ! real(dp), dimension(500*(Npol+1)):: dtrAng_s, dtxAng_s, dtyAng_s, rShape_s, rAng_s, xAng_s, yAng_s
    real(dp), dimension(500*(Npol+1)):: psi1_s, nu1_s, dchidx_s, dB_drho_s, dB_dl_s, dnu_drho_s, dnu_dl_s, grad_nu_s
    real(dp), dimension(500*(Npol+1)):: gxx, gxy, gxz, gyy, gyz, gzz, dtheta_dchi_s, dBp_dchi_s, jacobian, dBdx, dBdz
    real(dp), dimension(500*(Npol+1)):: g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, tmp_arr_s, dxdR_s, dxdZ_s, K_x, K_y !tmp_arr2

    real(dp), dimension(1:Nz):: gxx_out,gxy_out,gxz_out,gyy_out,gyz_out,gzz_out,Bfield_out,jacobian_out, dBdx_out, dBdz_out, chi_out
    real(dp), dimension(1:Nz):: R_out, Z_out, dxdR_out, dxdZ_out
    real(dp):: d_inv, drPsi, dxPsi, dq_dx, dq_dpsi, R0, Z0, B0, F, D0_full, D1_full, D2_full, D3_full
    !real(dp) :: Lnorm, Psi0 ! currently module-wide defined anyway
    real(dp):: pprime, ffprime, D0_mid, D1_mid, D2_mid, D3_mid, dx_drho, pi, mu_0, dzprimedz
    ! real(dp):: rho_a, psiN, drpsiN, CN2, CN3, Rcenter, Zcenter, drRcenter, drZcenter
    logical:: bMaxShift
    integer:: i, k, iBmax

    Npol_ext = Npol+2
    Npol_s   = Npol+1
    np   = 500*Npol_ext
    np_s = 500*Npol_s

    rho = trpeps*major_R
    if (rho.le.0.0) stop 'flux surface radius not defined'
    trpeps = rho/major_R

    q0 = sign_Ip_CW * sign_Bt_CW * abs(q0)

    R0=major_R
    B0=1.0*sign_Bt_CW
    F=R0*B0
    Z0=major_Z
    pi = acos(-1.0)
    mu_0=4.0*pi

    theta=linspace(-pi*Npol_ext,pi*Npol_ext-2._dp*pi*Npol_ext/np,np)
    d_inv=asin(delta)

    thetaShift = 0.0
    iBmax = 1

    !flux surface parametrization
    thAdj = theta + thetaShift
    if (zeta/=0.0 .or. s_zeta/=0.0) then
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

    !dl/dtheta
    dlp=(R_theta**2+Z_theta**2)**0.5

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
    drPsi=sign_Ip_CW*F/(2.*pi*Npol_ext*q0)*sum(tmp_arr)*2*pi*Npol_ext/np !dlp_int(tmp_arr,1.0)

    !Poloidal field (Bp = Bvec * nabla l)
    Bp=sign_Ip_CW * drPsi / J_r * dlp

    !toroidal field
    Bphi=F/R

    !total modulus of Bfield
    B=sqrt(Bphi**2+Bp**2)

    bMaxShift = .false.
    ! if (thetaShift==0.0.and.trim(magn_geometry).ne.'miller_general') then
    if (thetaShift==0.0) then
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

    !definition of radial coordinate! dx_drho=1 --> x = r
    dx_drho=1. !drPsi/Psi0*Lnorm*q0
    if (my_id==0) write(*,"(A,ES12.4)") 'Using radial coordinate with dx/dr = ',dx_drho
    dxPsi  = drPsi/dx_drho
    C_y    = dxPsi*sign_Ip_CW
    C_xy   = abs(B0*dxPsi/C_y)

    if (my_id==0) then
       write(*,"(A,ES12.4,A,ES12.4,A,ES12.4)") &
            "Setting C_xy = ",C_xy,' C_y = ', C_y," C_x' = ", 1./dxPsi
       write(*,'(A,ES12.4)') "B_unit/Bref conversion factor = ", q0/rho*drPsi
       write(*,'(A,ES12.4)') "dPsi/dr = ", drPsi
       if (thetaShift.ne.0.0) write(*,'(A,ES12.4)') "thetaShift = ", thetaShift
    endif


    !--------shear is expected to be defined as rho/q*dq/drho--------!
    dq_dx=shat*q0/rho/dx_drho
    dq_dpsi=dq_dx/dxPsi
    pprime=-amhd/q0**2/R0/(2*mu_0)*B0**2/drPsi

    !neg. dpdx normalized to magnetic pressure for pressure term
    dpdx_pm_geom=amhd/q0**2/R0/dx_drho

    !first coefficient of psi in varrho expansion
    psi1 = R*Bp*sign_Ip_CW

    !integrals for ffprime evaluation
    do i=1,np
       tmp_arr=(2./Rc-2.*cosu/R)/(R*psi1**2)
       D0(i)=-F*dlp_int_ind(tmp_arr,dlp,i)
       tmp_arr=B**2*R/psi1**3
       D1(i)=-dlp_int_ind(tmp_arr,dlp,i)/F
       tmp_arr=mu_0*R/psi1**3
       D2(i)=-dlp_int_ind(tmp_arr,dlp,i)*F
       tmp_arr=1./(R*psi1)
       D3(i)=-dlp_int_ind(tmp_arr,dlp,i)*F
    enddo
    tmp_arr=(2./Rc-2.*cosu/R)/(R*psi1**2)
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

    ffprime=-(sign_Ip_CW*dq_dpsi*2.*pi*Npol_ext+D0_full+D2_full*pprime)/D1_full

    if (my_id==0) then
       write(*,'(A,ES12.4)') "ffprime = ", ffprime
    endif
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

    if (sign_Ip_CW.lt.0.0) then !make chi increasing function to not confuse lag3interp
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

    if (zeta/=0.0 .or. s_zeta/=0.0) then
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
    if (sign_Ip_CW.lt.0.0) then
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

    !contravariant metric coefficients (varrho,l,fz)->(x,y,z)
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

    !curvature terms (these are just local and will be recalculated in geometry.F90)
    K_x = (0.-g_yz/g_zz*dBdz)
    K_y = (dBdx-g_xz/g_zz*dBdz)

    !(R,Z) derivatives for visualization
    dxdR_s = dx_drho/drPsi*psi1_s*cosu_s
    dxdZ_s = dx_drho/drPsi*psi1_s*sinu_s

    if (edge_opt==0.0) then
       !gene z-grid
       chi_out=linspace(-pi*Npol,pi*Npol-2*pi*Npol/Nz,Nz)
    else
       !new parallel coordinate chi_out==zprime
       !see also tracer_aux.F90
       if (Npol>1) STOP "ERROR: Npol>1 has not been implemented for edge_opt=\=0.0"
       do k=izs,ize
          chi_out(k)=sinh((-pi+k*2.*pi/Nz)*log(edge_opt*pi+sqrt(edge_opt**2*pi**2+1))/pi)/edge_opt
       enddo
       !transform metrics according to chain rule
       do k=1,np_s
         !>dz'/dz conversion for edge_opt as function of z
          if (edge_opt.gt.0) then
             dzprimedz = edge_opt*pi/log(edge_opt*pi+sqrt((edge_opt*pi)**2+1))/&
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

    !interpolate down to GENE z-grid
    call lag3interp(gxx,chi_s,np_s,gxx_out,chi_out,Nz)
    call lag3interp(gxy,chi_s,np_s,gxy_out,chi_out,Nz)
    call lag3interp(gxz,chi_s,np_s,gxz_out,chi_out,Nz)
    call lag3interp(gyy,chi_s,np_s,gyy_out,chi_out,Nz)
    call lag3interp(gyz,chi_s,np_s,gyz_out,chi_out,Nz)
    call lag3interp(gzz,chi_s,np_s,gzz_out,chi_out,Nz)
    call lag3interp(B_s,chi_s,np_s,Bfield_out,chi_out,Nz)
    call lag3interp(jacobian,chi_s,np_s,jacobian_out,chi_out,Nz)
    call lag3interp(dBdx,chi_s,np_s,dBdx_out,chi_out,Nz)
    call lag3interp(dBdz,chi_s,np_s,dBdz_out,chi_out,Nz)
    call lag3interp(R_s,chi_s,np_s,R_out,chi_out,Nz)
    call lag3interp(Z_s,chi_s,np_s,Z_out,chi_out,Nz)
    call lag3interp(dxdR_s,chi_s,np_s,dxdR_out,chi_out,Nz)
    call lag3interp(dxdZ_s,chi_s,np_s,dxdZ_out,chi_out,Nz)
    ! Fill the geom arrays with the results
    do eo=0,1
    gxx_(izs:ize,eo)      =gxx_out(izs:ize)
    gyy_(izs:ize,eo)      =gyy_out(izs:ize)
    gxz_(izs:ize,eo)      =gxz_out(izs:ize)
    gyz_(izs:ize,eo)      =gyz_out(izs:ize)
    dBdx_(izs:ize,eo)     =dBdx_out(izs:ize)
    dBdy_(izs:ize,eo)     =0.
    gxy_(izs:ize,eo)      =gxy_out(izs:ize)
    gzz_(izs:ize,eo)      =gzz_out(izs:ize)
    Bfield_(izs:ize,eo)   =Bfield_out(izs:ize)
    jacobian_(izs:ize,eo) =jacobian_out(izs:ize)
    dBdz_(izs:ize,eo)     =dBdz_out(izs:ize)
    R_hat_(izs:ize,eo)    =R_out(izs:ize)
    Z_hat_(izs:ize,eo)    =Z_out(izs:ize)
    dxdR_(izs:ize,eo)     = dxdR_out(izs:ize)
    dxdZ_(izs:ize,eo)     = dxdZ_out(izs:ize)

    !! UPDATE GHOSTS VALUES (since the miller function in GENE does not)
    CALL update_ghosts_z(gxx_(:,eo))
    CALL update_ghosts_z(gyy_(:,eo))
    CALL update_ghosts_z(gxz_(:,eo))
    CALL update_ghosts_z(gxy_(:,eo))
    CALL update_ghosts_z(gzz_(:,eo))
    CALL update_ghosts_z(Bfield_(:,eo))
    CALL update_ghosts_z(dBdx_(:,eo))
    CALL update_ghosts_z(dBdy_(:,eo))
    CALL update_ghosts_z(dBdz_(:,eo))
    CALL update_ghosts_z(jacobian_(:,eo))
    CALL update_ghosts_z(R_hat_(:,eo))
    CALL update_ghosts_z(Z_hat_(:,eo))
    CALL update_ghosts_z(dxdR_(:,eo))
    CALL update_ghosts_z(dxdZ_(:,eo))

    ! Curvature operator (Frei et al. 2022 eq 2.15)
    DO iz = izgs,izge
      G1 = gxy_(iz,eo)*gxy_(iz,eo)-gxx_(iz,eo)*gyy_(iz,eo)
      G2 = gxy_(iz,eo)*gxz_(iz,eo)-gxx_(iz,eo)*gyz_(iz,eo)
      G3 = gyy_(iz,eo)*gxz_(iz,eo)-gxy_(iz,eo)*gyz_(iz,eo)
      Cx = (G1*dBdy_(iz,eo) + G2*dBdz_(iz,eo))/Bfield_(iz,eo)
      Cy = (G3*dBdz_(iz,eo) - G1*dBdx_(iz,eo))/Bfield_(iz,eo)

      DO iky = ikys, ikye
        ky = kyarray(iky)
         DO ikx= ikxs, ikxe
           kx = kxarray(ikx)
           Ckxky_(iky, ikx, iz,eo) = (Cx*kx + Cy*ky)
         ENDDO
      ENDDO
      ! coefficient in the front of parallel derivative
      gradz_coeff_(iz,eo) = 1._dp / jacobian_(iz,eo) / Bfield_(iz,eo)
    ENDDO
  ENDDO

  contains


    SUBROUTINE update_ghosts_z(fz_)
      IMPLICIT NONE
      ! INTEGER,  INTENT(IN) :: nztot_
      REAL(dp), DIMENSION(izgs:izge), INTENT(INOUT) :: fz_
      REAL(dp), DIMENSION(-2:2) :: buff
      INTEGER :: status(MPI_STATUS_SIZE), count

      IF(Nz .GT. 1) THEN
        IF (num_procs_z .GT. 1) THEN
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            count = 1 ! one point to exchange
            !!!!!!!!!!! Send ghost to up neighbour !!!!!!!!!!!!!!!!!!!!!!
            CALL mpi_sendrecv(fz_(ize), count, MPI_DOUBLE, nbr_U, 0, & ! Send to Up the last
                              buff(-1), count, MPI_DOUBLE, nbr_D, 0, & ! Receive from Down the first-1
                              comm0, status, ierr)

            CALL mpi_sendrecv(fz_(ize-1), count, MPI_DOUBLE, nbr_U, 0, & ! Send to Up the last
                                buff(-2), count, MPI_DOUBLE, nbr_D, 0, & ! Receive from Down the first-2
                              comm0, status, ierr)

            !!!!!!!!!!! Send ghost to down neighbour !!!!!!!!!!!!!!!!!!!!!!
            CALL mpi_sendrecv(fz_(izs), count, MPI_DOUBLE, nbr_D, 0, & ! Send to Down the first
                              buff(+1), count, MPI_DOUBLE, nbr_U, 0, & ! Recieve from Up the last+1
                              comm0, status, ierr)

            CALL mpi_sendrecv(fz_(izs+1), count, MPI_DOUBLE, nbr_D, 0, & ! Send to Down the first
                                buff(+2), count, MPI_DOUBLE, nbr_U, 0, & ! Recieve from Up the last+2
                              comm0, status, ierr)
         ELSE
           buff(-1) = fz_(ize  )
           buff(-2) = fz_(ize-1)
           buff(+1) = fz_(izs  )
           buff(+2) = fz_(izs+1)
         ENDIF
         fz_(ize+1) = buff(+1)
         fz_(ize+2) = buff(+2)
         fz_(izs-1) = buff(-1)
         fz_(izs-2) = buff(-2)
      ENDIF
    END SUBROUTINE update_ghosts_z


    !> Generate an equidistant array from min to max with n points
    function linspace(min,max,n) result(out)
      real(dp), INTENT(IN):: min, max
      integer,  INTENT(IN):: n
      real(dp), dimension(n):: out

      do i=1,n
         out(i)=min+(i-1)*(max-min)/(n-1)
      enddo
    end function linspace

    !> Weighted average
    real(dp) function average(var,weight)
      real(dp), dimension(np), INTENT(IN):: var, weight
      average=sum(var*weight)/sum(weight)
    end function average

    !> full theta integral with weight function dlp
    real(dp) function dlp_int(var,dlp)
      real(dp), dimension(np), INTENT(IN):: var, dlp
      dlp_int=sum(var*dlp)*2*pi*Npol_ext/np
    end function dlp_int

    !> theta integral with weight function dlp, up to index 'ind'
    real(dp) function dlp_int_ind(var,dlp,ind)
      real(dp), dimension(np), INTENT(IN):: var, dlp
      integer, INTENT(IN):: ind

      dlp_int_ind=0.
      if (ind.gt.1) then
         dlp_int_ind=dlp_int_ind+var(1)*dlp(1)*pi*Npol_ext/np
         dlp_int_ind=dlp_int_ind+(sum(var(2:ind-1)*dlp(2:ind-1)))*2*pi*Npol_ext/np
         dlp_int_ind=dlp_int_ind+var(ind)*dlp(ind)*pi*Npol_ext/np
      endif
    end function dlp_int_ind

    !> 1st derivative with 2nd order finite differences
    function deriv_fd(y,x,n) result(out)
      integer,                INTENT(IN) :: n
      real(dp), dimension(n), INTENT(IN):: x,y
      real(dp), dimension(n) :: out,dx

      !call lag3deriv(y,x,n,out,x,n)

      out=0.
      do i=2,n-1
         out(i)=out(i)-y(i-1)/2
         out(i)=out(i)+y(i+1)/2
      enddo
      out(1)=y(2)-y(1)
      out(n)=y(n)-y(n-1)
      dx=x(2)-x(1)
      out=out/dx

    end function deriv_fd


  end subroutine get_miller


END MODULE miller
