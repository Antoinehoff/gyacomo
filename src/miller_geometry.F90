#include "redef.h"
!! This source has been adapted from GENE https://genecode.org/ !!
!>Implementation of the local equilibrium model of [R.L. Miller et al., PoP 5, 973 (1998)
!>and [J. Candy, PPCF 51, 105009 (2009)]
MODULE miller_mod
  ! use coordinates,only: gcoor, get_dzprimedz
  ! use discretization
  use lagrange_interpolation
  ! use par_geom
  ! use par_in, only: beta, sign_Ip_CW, sign_Bt_CW, n_pol
  ! use par_other, only: print_ini_msg

  implicit none
  public:: get_miller, set_miller_defaults
  public:: rho, kappa, delta, s_kappa, s_delta, drR, drZ, zeta, s_zeta
  public:: thetaShift
  public:: mMode, nMode
  public:: thetak, thetad
  public:: aSurf, Delta2, Delta3, theta2, theta3, Raxis, Zaxis
  public:: Deltam, Deltan, s_Deltam, s_Deltan, thetam, thetan
  public:: cN_m, sN_m, cNdr_m, sNdr_m

  private

  REAL :: rho, kappa, delta, s_kappa, s_delta, drR, drZ, zeta, s_zeta

  INTEGER :: mMode, nMode
  REAL :: thetaShift
  REAL :: thetak, thetad
  REAL :: aSurf, Delta2, Delta3, theta2, theta3, Raxis, Zaxis
  REAL :: Deltam, Deltan, s_Deltam, s_Deltan, thetam, thetan

  INTEGER, PARAMETER :: IND_M=32
  REAL, DIMENSION(0:IND_M-1) :: cN_m, sN_m, cNdr_m, sNdr_m

CONTAINS

  !>Set defaults for miller parameters
  subroutine set_miller_defaults
    rho = -1.0
    kappa = 1.0
    s_kappa = 0.0
    delta = 0.0
    s_delta = 0.0
    drR = 0.0
    drZ = 0.0
    zeta = 0.0
    s_zeta = 0.0


    thetak = 0.0
    thetad = 0.0

    aSurf = 0.54
    Delta2 = 1.0
    Delta3 = 1.0
    theta2 = 0.0
    theta3 = 0.0
    Raxis = 1.0
    Zaxis = 0.0

    mMode = 2
    nMode = 3
    Deltam = 1.0
    Deltan = 1.0
    s_Deltam = 0.0
    s_Deltan = 0.0
    thetam = 0.0
    thetan = 0.0

    cN_m = 0.0
    sN_m = 0.0
    cNdr_m = 0.0
    sNdr_m = 0.0
  end subroutine set_miller_defaults

  !>Get Miller metric, magnetic field, jacobian etc.
  subroutine get_miller(geom,edge_opt)  !previously: get_miller_a
    type(geomtype),intent(inout):: geom
    real, intent(in):: edge_opt
    integer:: np, np_s, n_pol_ext, n_pol_s

    real, dimension(500*(n_pol+2)):: R,Z,R_rho,Z_rho,R_theta,Z_theta,R_theta_theta,Z_theta_theta,dlp,Rc,cosu,sinu,Bphi
    real, dimension(500*(n_pol+2)):: drRcirc, drRelong, drRelongTilt, drRtri, drRtriTilt, drZcirc, drZelong, drZelongTilt
    real, dimension(500*(n_pol+2)):: drZtri, drZtriTilt, dtdtRcirc, dtdtRelong, dtdtRelongTilt, dtdtRtri, dtdtRtriTilt
    real, dimension(500*(n_pol+2)):: dtdtZcirc, dtdtZelong, dtdtZelongTilt, dtdtZtri, dtdtZtriTilt, dtRcirc, dtRelong
    real, dimension(500*(n_pol+2)):: dtRelongTilt, dtRtri, dtRtriTilt, dtZcirc, dtZelong, dtZelongTilt, dtZtri, dtZtriTilt
    real, dimension(500*(n_pol+2)):: Rcirc, Relong, RelongTilt, Rtri, RtriTilt, Zcirc, Zelong, ZelongTilt, Ztri, ZtriTilt
    real, dimension(500*(n_pol+2)):: drrShape, drrAng, drxAng, dryAng, dtdtrShape, dtdtrAng, dtdtxAng
    real, dimension(500*(n_pol+2)):: dtdtyAng, dtrShape, dtrAng, dtxAng, dtyAng, rShape, rAng, xAng, yAng
    real, dimension(500*(n_pol+2)):: theta, thAdj, J_r, B, Bp, D0, D1, D2, D3, nu, chi, psi1, nu1
    real, dimension(500*(n_pol+2)):: tmp_reverse, theta_reverse, tmp_arr

    real, dimension(500*(n_pol+1)):: theta_s, thAdj_s, chi_s, theta_s_reverse
    real, dimension(500*(n_pol+1)):: R_s, Z_s, R_theta_s, Z_theta_s, Rc_s, cosu_s, sinu_s, Bphi_s, B_s, Bp_s, dlp_s
    real, dimension(500*(n_pol+1)):: dtRcirc_s, dtRelong_s, dtRelongTilt_s, dtRtri_s, dtRtriTilt_s, dtZcirc_s
    real, dimension(500*(n_pol+1)):: dtZelong_s, dtZelongTilt_s, dtZtri_s, dtZtriTilt_s, Rcirc_s, Relong_s, RelongTilt_s
    real, dimension(500*(n_pol+1)):: Rtri_s, RtriTilt_s, Zcirc_s, Zelong_s, ZelongTilt_s, Ztri_s, ZtriTilt_s, dtrShape_s
    real, dimension(500*(n_pol+1)):: dtrAng_s, dtxAng_s, dtyAng_s, rShape_s, rAng_s, xAng_s, yAng_s
    real, dimension(500*(n_pol+1)):: psi1_s, nu1_s, dchidx_s, dB_drho_s, dB_dl_s, dnu_drho_s, dnu_dl_s, grad_nu_s
    real, dimension(500*(n_pol+1)):: gxx, gxy, gxz, gyy, gyz, gzz, dtheta_dchi_s, dBp_dchi_s, jacobian, dBdx, dBdz
    real, dimension(500*(n_pol+1)):: g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, tmp_arr_s, dxdR_s, dxdZ_s, K_x, K_y !tmp_arr2

    real, dimension(0:nz0-1):: gxx_out,gxy_out,gxz_out,gyy_out,gyz_out,gzz_out,Bfield_out,jacobian_out, dBdx_out, dBdz_out, chi_out
    real, dimension(0:nz0-1):: R_out, Z_out, dxdR_out, dxdZ_out
    real:: d_inv, drPsi, dxPsi, dq_dx, dq_dpsi, R0, Z0, B0, F, D0_full, D1_full, D2_full, D3_full
    !real :: Lnorm, Psi0 ! currently module-wide defined anyway
    real:: pprime, ffprime, D0_mid, D1_mid, D2_mid, D3_mid, dx_drho, pi, mu_0, dzprimedz
    real:: rho_a, psiN, drpsiN, CN2, CN3, Rcenter, Zcenter, drRcenter, drZcenter
    logical:: bMaxShift
    integer:: i, k, iBmax

    n_pol_ext = n_pol+2
    n_pol_s   = n_pol+1
    np   = 500*n_pol_ext
    np_s = 500*n_pol_s

    if (rho.lt.0.0) rho = trpeps*major_R
    if (rho.le.0.0) stop 'flux surface radius not defined'
    trpeps = rho/major_R
    gcoor%x0=rho

    q0 = sign_Ip_CW * sign_Bt_CW * abs(q0)

    R0=major_R
    B0=1.0*sign_Bt_CW
    F=R0*B0
    Z0=major_Z
    pi = acos(-1.0)
    mu_0=4.0*pi

    theta=linspace(-pi*n_pol_ext,pi*n_pol_ext-2*pi*n_pol_ext/np,np)
    d_inv=asin(delta)

    thetaShift = 0.0
    iBmax = 1

    !flux surface parametrization
    thAdj = theta + thetaShift
    select case (magn_geometry)
    case ('miller')
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
    case default
      write (*,*) "ERROR: invalid analytic geometry specification"
    end select

    !dl/dtheta
    dlp=(R_theta**2+Z_theta**2)**0.5

    !curvature radius
    Rc=dlp**3*(R_theta*Z_theta_theta-Z_theta*R_theta_theta)**(-1)

    ! some useful quantities (see papers for definition of u)
    cosu=Z_theta/dlp
    sinu=-R_theta/dlp

    !Jacobian J_r = (dPsi/dr) J_psi = (dPsi/dr) / [(nabla phi x nabla psi)* nabla theta]
    !             = R * (dR/drho dZ/dtheta - dR/dtheta dZ/drho) = R dlp / |nabla r|
    J_r=R*(R_rho*Z_theta-R_theta*Z_rho)

    !From definition of q = 1/(2 pi) int (B nabla phi) / (B nabla theta) dtheta:
    !dPsi/dr = sign_Bt sign_Ip / (2 pi q) int F / R^2 J_r dtheta
    !        = F / (2 pi |q|) int J_r/R^2 dtheta
    tmp_arr=J_r/R**2
    drPsi=sign_Ip_CW*F/(2.*pi*n_pol_ext*q0)*sum(tmp_arr)*2*pi*n_pol_ext/np !dlp_int(tmp_arr,1.0)

    !Poloidal field (Bp = Bvec * nabla l)
    Bp=sign_Ip_CW * drPsi / J_r * dlp

    !toroidal field
    Bphi=F/R

    !total modulus of Bfield
    B=sqrt(Bphi**2+Bp**2)

    bMaxShift = .false.
    if (thetaShift==0.0.and.trim(magn_geometry).ne.'miller_general') then
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
    if (mype==0.and.print_ini_msg) write(*,"(A,ES12.4)") 'Using radial coordinate with dx/dr = ',dx_drho
    dxPsi=drPsi/dx_drho
    geom%C_y=dxPsi*sign_Ip_CW
    geom%Cyq0_x0=geom%C_y(gpdisc%pi1gl)*q0/rho
    geom%C_xy=abs(B0*dxPsi/geom%C_y)

    if (mype==0.and.print_ini_msg) then
       write(*,"(A,ES12.4,A,ES12.4,A,ES12.4)") &
            "Setting C_xy = ",geom%C_xy(gpdisc%pi1gl),' C_y = ', geom%C_y(gpdisc%pi1gl)," C_x' = ", 1./dxPsi
       write(*,'(A,ES12.4)') "B_unit/Bref conversion factor = ", q0/rho*drPsi
       write(*,'(A,ES12.4)') "dPsi/dr = ", drPsi
       if (thetaShift.ne.0.0) write(*,'(A,ES12.4)') "thetaShift = ", thetaShift
    endif


    !--------shear is expected to be defined as rho/q*dq/drho--------!
    dq_dx=shat*q0/rho/dx_drho
    dq_dpsi=dq_dx/dxPsi
    pprime=-amhd/q0**2/R0/(2*mu_0)*B0**2/drPsi

    !neg. dpdx normalized to magnetic pressure for pressure term
    geom%dpdx_pm_geom=amhd/q0**2/R0/dx_drho

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

    ffprime=-(sign_Ip_CW*dq_dpsi*2.*pi*n_pol_ext+D0_full+D2_full*pprime)/D1_full

    if (mype==0.and.print_ini_msg) then
       write(*,'(A,ES12.4)') "ffprime = ", ffprime
    endif
    D0=D0-D0_mid
    D1=D1-D1_mid
    D2=D2-D2_mid
    nu=D3-D3_mid

    nu1=psi1*(D0+D1*ffprime+D2*pprime)

    !straight field line angle defined on equidistant theta grid
    !alpha = phi + nu = - (q chi - phi) => chi = -nu / q
    chi=-nu/q0

    !correct small scaling error (<0.5%, due to finite integration resolution)
    chi=chi*(maxval(theta)-minval(theta))/(maxval(chi)-minval(chi))

    !new grid equidistant in straight field line angle
    chi_s = linspace(-pi*n_pol_s,pi*n_pol_s-2*pi*n_pol_s/np_s,np_s)

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

    select case (magn_geometry)
    case ('miller')
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
    case default
       write (*,*) "ERROR: invalid analytic geometry specification"
    end select

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

    !Bfield derivatives in Mercier-Luc coordinates (varrho,l,phi)
    dB_drho_s=-1./B_s*(F**2/R_s**3*cosu_s+Bp_s**2/Rc_s+mu_0*psi1_s*pprime)
    dB_dl_s=1./B_s*(Bp_s*dBp_dchi_s/dtheta_dchi_s/dlp_s+F**2/R_s**3*sinu_s)

    dnu_drho_s=nu1_s
    dnu_dl_s=-F/(R_s*psi1_s)
    grad_nu_s=sqrt(dnu_drho_s**2+dnu_dl_s**2)

    !contravariant metric coefficients (varrho,l,phi)->(x,y,z)
    gxx=(psi1_s/dxPsi)**2
    gxy=-psi1_s/dxPsi*geom%C_y(gpdisc%pi1gl)*sign_Ip_CW*nu1_s
    gxz=-psi1_s/dxPsi*(nu1_s+psi1_s*dq_dpsi*chi_s)/q0
    gyy=geom%C_y(gpdisc%pi1gl)**2*(grad_nu_s**2+1/R_s**2)
    gyz=sign_Ip_CW*geom%C_y(gpdisc%pi1gl)/q0*(grad_nu_s**2+dq_dpsi*nu1_s*psi1_s*chi_s)
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
    dBdx=jacobian*geom%C_y(gpdisc%pi1gl)/(q0*R_s)*(F/(R_s*psi1_s)*dB_drho_s+(nu1_s+dq_dpsi*chi_s*psi1_s)*dB_dl_s)
    dBdz=1./B_s*(Bp_s*dBp_dchi_s-F**2/R_s**3*R_theta_s)

    !curvature terms (these are just local and will be recalculated in geometry.F90)
    K_x = (0.-g_yz/g_zz*dBdz)
    K_y = (dBdx-g_xz/g_zz*dBdz)

    !(R,Z) derivatives for visualization
    dxdR_s = dx_drho/drPsi*psi1_s*cosu_s
    dxdZ_s = dx_drho/drPsi*psi1_s*sinu_s

    if (edge_opt==0.0) then
       !gene z-grid
       chi_out=linspace(-pi*n_pol,pi*n_pol-2*pi*n_pol/nz0,nz0)
    else
       !new parallel coordinate chi_out==zprime
       !see also tracer_aux.F90
       if (n_pol>1) STOP "ERROR: n_pol>1 has not been implemented for edge_opt=\=0.0"
       do k=0,nz0-1
          chi_out(k)=sinh((-pi+k*2.*pi/nz0)*log(edge_opt*pi+sqrt(edge_opt**2*pi**2+1))/pi)/edge_opt
       enddo
       !transform metrics according to chain rule
       do k=1,np_s
          dzprimedz=get_dzprimedz(chi_s(k))
          gxz(k)=gxz(k)*dzprimedz
          gyz(k)=gyz(k)*dzprimedz
          gzz(k)=gzz(k)*dzprimedz**2
          jacobian(k)=jacobian(k)/dzprimedz
          dBdz(k)=dBdz(k)/dzprimedz
       enddo
    endif !edge_opt

    !interpolate down to GENE z-grid
    call lag3interp(gxx,chi_s,np_s,gxx_out,chi_out,nz0)
    call lag3interp(gxy,chi_s,np_s,gxy_out,chi_out,nz0)
    call lag3interp(gxz,chi_s,np_s,gxz_out,chi_out,nz0)
    call lag3interp(gyy,chi_s,np_s,gyy_out,chi_out,nz0)
    call lag3interp(gyz,chi_s,np_s,gyz_out,chi_out,nz0)
    call lag3interp(gzz,chi_s,np_s,gzz_out,chi_out,nz0)
    call lag3interp(B_s,chi_s,np_s,Bfield_out,chi_out,nz0)
    call lag3interp(jacobian,chi_s,np_s,jacobian_out,chi_out,nz0)
    call lag3interp(dBdx,chi_s,np_s,dBdx_out,chi_out,nz0)
    call lag3interp(dBdz,chi_s,np_s,dBdz_out,chi_out,nz0)
    call lag3interp(R_s,chi_s,np_s,R_out,chi_out,nz0)
    call lag3interp(Z_s,chi_s,np_s,Z_out,chi_out,nz0)
    call lag3interp(dxdR_s,chi_s,np_s,dxdR_out,chi_out,nz0)
    call lag3interp(dxdZ_s,chi_s,np_s,dxdZ_out,chi_out,nz0)

    !select local k range
    do i=gpdisc%pi1gl,gpdisc%pi2gl
       if (gdisc%yx_order) then
          geom%gii(i,lk1:lk2)=gyy_out(lk1:lk2)
          geom%gjj(i,lk1:lk2)=gxx_out(lk1:lk2)
          geom%giz(i,lk1:lk2)=gyz_out(lk1:lk2)
          geom%gjz(i,lk1:lk2)=gxz_out(lk1:lk2)
          geom%dBdi(i,lk1:lk2)=0.
          geom%dBdj(i,lk1:lk2)=dBdx_out(lk1:lk2)
       else
          geom%gii(i,lk1:lk2)=gxx_out(lk1:lk2)
          geom%gjj(i,lk1:lk2)=gyy_out(lk1:lk2)
          geom%giz(i,lk1:lk2)=gxz_out(lk1:lk2)
          geom%gjz(i,lk1:lk2)=gyz_out(lk1:lk2)
          geom%dBdi(i,lk1:lk2)=dBdx_out(lk1:lk2)
          geom%dBdj(i,lk1:lk2)=0.
       endif
       geom%gij(i,lk1:lk2)=gxy_out(lk1:lk2)
       geom%gzz(i,lk1:lk2)=gzz_out(lk1:lk2)
       geom%Bfield(i,lk1:lk2)=Bfield_out(lk1:lk2)
       geom%jacobian(i,lk1:lk2)=jacobian_out(lk1:lk2)
       geom%dBdz(i,lk1:lk2)=dBdz_out(lk1:lk2)
       if (Lref.gt.0.) then
          geom%R(i,lk1:lk2)=R_out(lk1:lk2)*Lref
          geom%Z(i,lk1:lk2)=Z_out(lk1:lk2)*Lref
       else
          geom%R(i,lk1:lk2)=R_out(lk1:lk2)
          geom%Z(i,lk1:lk2)=Z_out(lk1:lk2)
       endif
       geom%R_hat(i,lk1:lk2)=R_out(lk1:lk2)
       geom%Z_hat(i,lk1:lk2)=Z_out(lk1:lk2)
       geom%dxdR(i,lk1:lk2)= dxdR_out(lk1:lk2)
       geom%dxdZ(i,lk1:lk2)= dxdZ_out(lk1:lk2)
    enddo

    geom%q_prof=q0
    geom%dqdx_prof=shat*q0/gcoor%x0

  contains

    !> Generate an equidistant array from min to max with n points
    function linspace(min,max,n) result(out)
      real:: min, max
      integer:: n
      real, dimension(n):: out

      do i=1,n
         out(i)=min+(i-1)*(max-min)/(n-1)
      enddo
    end function linspace

    !> Weighted average
    real function average(var,weight)
      real, dimension(np):: var, weight
      average=sum(var*weight)/sum(weight)
    end function average

    !> full theta integral with weight function dlp
    real function dlp_int(var,dlp)
      real, dimension(np):: var, dlp
      dlp_int=sum(var*dlp)*2*pi*n_pol_ext/np
    end function dlp_int

    !> theta integral with weight function dlp, up to index 'ind'
    real function dlp_int_ind(var,dlp,ind)
      real, dimension(np):: var, dlp
      integer:: ind

      dlp_int_ind=0.
      if (ind.gt.1) then
         dlp_int_ind=dlp_int_ind+var(1)*dlp(1)*pi*n_pol_ext/np
         dlp_int_ind=dlp_int_ind+(sum(var(2:ind-1)*dlp(2:ind-1)))*2*pi*n_pol_ext/np
         dlp_int_ind=dlp_int_ind+var(ind)*dlp(ind)*pi*n_pol_ext/np
      endif
    end function dlp_int_ind

    !> 1st derivative with 2nd order finite differences
    function deriv_fd(y,x,n) result(out)
      integer, intent(in) :: n
      real, dimension(n):: x,y,out,dx

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


END MODULE miller_mod
