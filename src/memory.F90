SUBROUTINE memory
  ! Allocate arrays (done dynamically otherwise size is unknown)

  USE array
  USE basic, ONLY: allocate_array
  USE fields
  USE grid, ONLY: local_Na, local_Np,Ngp ,local_Nj,Ngj, local_nky, local_nkx,local_Nz,Ngz, jmax, pmax
  USE collision
  USE time_integration, ONLY: ntimelevel
  USE prec_const
  USE model, ONLY: Na
  IMPLICIT NONE

  ! Electrostatic potential
  CALL allocate_array(           phi, 1,local_nky, 1,local_nkx, 1,local_Nz+Ngz)
  CALL allocate_array(           psi, 1,local_nky, 1,local_nkx, 1,local_Nz+Ngz)
  CALL allocate_array(inv_poisson_op, 1,local_nky, 1,local_nkx, 1,local_Nz)
  CALL allocate_array( inv_ampere_op, 1,local_nky, 1,local_nkx, 1,local_Nz)
  CALL allocate_array(   inv_pol_ion, 1,local_nky, 1,local_nkx, 1,local_Nz)
  ! CALL allocate_array(HF_phi_correction_operator, 1,local_nky, 1,local_nkx, 1,local_Nz)

  !Moments related arrays
  CALL allocate_array(           Na00, 1,local_Na, 1,local_nky, 1,local_nkx, 1,local_Nz)
  CALL allocate_array(           dens, 1,local_Na, 1,local_nky, 1,local_nkx, 1,local_Nz)
  CALL allocate_array(           upar, 1,local_Na, 1,local_nky, 1,local_nkx, 1,local_Nz)
  CALL allocate_array(           uper, 1,local_Na, 1,local_nky, 1,local_nkx, 1,local_Nz)
  CALL allocate_array(           Tpar, 1,local_Na, 1,local_nky, 1,local_nkx, 1,local_Nz)
  CALL allocate_array(           Tper, 1,local_Na, 1,local_nky, 1,local_nkx, 1,local_Nz)
  CALL allocate_array(           temp, 1,local_Na, 1,local_nky, 1,local_nkx, 1,local_Nz)
  CALL allocate_array(         Kernel, 1,local_Na,                 1,local_Nj+Ngj, 1,local_nky, 1,local_nkx, 1,local_Nz+Ngz, 1,2)
  CALL allocate_array(        moments, 1,local_Na, 1,local_Np+Ngp, 1,local_Nj+Ngj, 1,local_nky, 1,local_nkx, 1,local_Nz+Ngz, 1,ntimelevel )
  CALL allocate_array(          Napjz, 1,local_Na, 1,local_Np,     1,local_Nj,                               1,local_Nz)
  CALL allocate_array(    moments_rhs, 1,local_Na, 1,local_Np,     1,local_Nj,     1,local_nky, 1,local_nkx, 1,local_Nz,     1,ntimelevel )
  CALL allocate_array( nadiab_moments, 1,local_Na, 1,local_Np+Ngp, 1,local_Nj+Ngj, 1,local_nky, 1,local_nkx, 1,local_Nz+Ngz)
  CALL allocate_array(       ddz_napj, 1,local_Na, 1,local_Np+Ngp, 1,local_Nj+Ngj, 1,local_nky, 1,local_nkx, 1,local_Nz)
  CALL allocate_array(     ddzND_Napj, 1,local_Na, 1,local_Np+Ngp, 1,local_Nj+Ngj, 1,local_nky, 1,local_nkx, 1,local_Nz)
  CALL allocate_array(    interp_napj, 1,local_Na, 1,local_Np+Ngp, 1,local_Nj+Ngj, 1,local_nky, 1,local_nkx, 1,local_Nz)
  CALL allocate_array(          Capj,  1,local_Na, 1,local_Np,     1,local_Nj,     1,local_nky, 1,local_nkx, 1,local_Nz)
  CALL allocate_array(          Sapj,  1,local_Na, 1,local_Np,     1,local_Nj,     1,local_nky, 1,local_nkx, 1,local_Nz)
  CALL allocate_array(     xnapj, 1,local_Na, 1,local_Np, 1,local_Nj)
  CALL allocate_array(   xnapp2j, 1,local_Na, 1,local_Np)
  CALL allocate_array(   xnapp1j, 1,local_Na, 1,local_Np)
  CALL allocate_array(   xnapm1j, 1,local_Na, 1,local_Np)
  CALL allocate_array(   xnapm2j, 1,local_Na, 1,local_Np)
  CALL allocate_array(   xnapjp1, 1,local_Na, 1,local_Nj)
  CALL allocate_array(   xnapjm1, 1,local_Na, 1,local_Nj)
  CALL allocate_array(   ynapp1j, 1,local_Na, 1,local_Np, 1,local_Nj)
  CALL allocate_array(   ynapm1j, 1,local_Na, 1,local_Np, 1,local_Nj)
  CALL allocate_array( ynapp1jm1, 1,local_Na, 1,local_Np, 1,local_Nj)
  CALL allocate_array( ynapm1jm1, 1,local_Na, 1,local_Np, 1,local_Nj)
  CALL allocate_array(   zNapm1j, 1,local_Na, 1,local_Np, 1,local_Nj)
  CALL allocate_array( zNapm1jp1, 1,local_Na, 1,local_Np, 1,local_Nj)
  CALL allocate_array( zNapm1jm1, 1,local_Na, 1,local_Np, 1,local_Nj)

  ! Non linear terms and dnjs table
  CALL allocate_array( dnjs, 1,jmax+1, 1,jmax+1, 1,jmax+1)

  ! Hermite fourth derivative coeff table 4*sqrt(p!/(p-4)!)
  CALL allocate_array( dv4_Hp_coeff, -2, pmax)

  CALL allocate_array( xphij,   1,local_Na, 1,local_Np, 1,local_Nj)
  CALL allocate_array( xphijp1, 1,local_Na, 1,local_Np, 1,local_Nj)
  CALL allocate_array( xphijm1, 1,local_Na, 1,local_Np, 1,local_Nj)
  CALL allocate_array( xpsij,   1,local_Na, 1,local_Np, 1,local_Nj)
  CALL allocate_array( xpsijp1, 1,local_Na, 1,local_Np, 1,local_Nj)
  CALL allocate_array( xpsijm1, 1,local_Na, 1,local_Np, 1,local_Nj)

  !___________________ 2x5D ARRAYS __________________________
  !! Collision matrices
  IF (GK_CO) THEN !GK collision matrices (one for each kperp)
    CALL allocate_array(  Cab_F, 1,Na, 1,Na, 1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), 1,local_nky, 1,local_nkx, 1,local_Nz)
    CALL allocate_array(  Cab_T, 1,Na, 1,Na, 1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), 1,local_nky, 1,local_nkx, 1,local_Nz)
    CALL allocate_array(  Caa,   1,Na,       1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), 1,local_nky, 1,local_nkx, 1,local_Nz)
  ELSE !DK collision matrix (same for every k)
      CALL allocate_array(  Cab_F, 1,Na, 1,Na, 1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), 1,1, 1,1, 1,1)
      CALL allocate_array(  Cab_T, 1,Na, 1,Na, 1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), 1,1, 1,1, 1,1)
      CALL allocate_array(  Caa,    1,Na,      1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), 1,1, 1,1, 1,1)
ENDIF

END SUBROUTINE memory
