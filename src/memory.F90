SUBROUTINE memory
  ! Allocate arrays (done dynamically otherwise size is unknown)

  USE array
  USE basic, ONLY: allocate_array
  USE fields
  USE grid, ONLY: local_na, local_np,ngp ,local_nj,ngj, local_nky, local_nkx,local_nz,ngz, jmax, pmax, nzgrid
  USE collision
  USE time_integration, ONLY: ntimelevel
  USE prec_const
  USE model, ONLY: na
  IMPLICIT NONE
  !==================== Ghosted arrays ====================
  ! 5D+ arrays
  CALL allocate_array(        moments, 1,local_na, 1,local_np+ngp, 1,local_nj+ngj, 1,local_nky, 1,local_nkx, 1,local_nz+ngz, 1,ntimelevel )
  CALL allocate_array( nadiab_moments, 1,local_na, 1,local_np+ngp, 1,local_nj+ngj, 1,local_nky, 1,local_nkx, 1,local_nz+ngz)
  CALL allocate_array(       ddz_napj, 1,local_na, 1,local_np+ngp, 1,local_nj+ngj, 1,local_nky, 1,local_nkx, 1,local_nz)
  CALL allocate_array(     ddzND_napj, 1,local_na, 1,local_np+ngp, 1,local_nj+ngj, 1,local_nky, 1,local_nkx, 1,local_nz)
  CALL allocate_array(    interp_napj, 1,local_na, 1,local_np+ngp, 1,local_nj+ngj, 1,local_nky, 1,local_nkx, 1,local_nz)
  
  ! 4D+ arrays
  CALL allocate_array(     Kernel, 1,local_na, 1,local_nj+ngj, 1,local_nky, 1,local_nkx, 1,local_nz+ngz, 1,nzgrid)
  
  ! 3D arrays
  CALL allocate_array( phi, 1,local_nky, 1,local_nkx, 1,local_nz+ngz)
  CALL allocate_array( psi, 1,local_nky, 1,local_nkx, 1,local_nz+ngz)
  
  ! smaller arrays
  
  
  !==================== Not ghosted arrays ====================
  ! 5D+ arrays
  CALL allocate_array(   moments_rhs, 1,local_na, 1,local_np,  1,local_nj, 1,local_nky, 1,local_nkx, 1,local_nz, 1,ntimelevel )
  CALL allocate_array(          Capj, 1,local_na, 1,local_np,  1,local_nj, 1,local_nky, 1,local_nkx, 1,local_nz)
  CALL allocate_array(          Sapj, 1,local_na, 1,local_np,  1,local_nj, 1,local_nky, 1,local_nkx, 1,local_nz)
  
  ! 4D+ arrays
  CALL allocate_array(         napjz, 1,local_na,  1,local_np,  1,local_nj, 1,local_nz)
  CALL allocate_array(          dens, 1,local_na, 1,local_nky, 1,local_nkx, 1,local_nz)
  CALL allocate_array(          upar, 1,local_na, 1,local_nky, 1,local_nkx, 1,local_nz)
  CALL allocate_array(          uper, 1,local_na, 1,local_nky, 1,local_nkx, 1,local_nz)
  CALL allocate_array(          Tpar, 1,local_na, 1,local_nky, 1,local_nkx, 1,local_nz)
  CALL allocate_array(          Tper, 1,local_na, 1,local_nky, 1,local_nkx, 1,local_nz)
  CALL allocate_array(          temp, 1,local_na, 1,local_nky, 1,local_nkx, 1,local_nz)
  
  ! 3D arrays
  CALL allocate_array(inv_poisson_op, 1,local_nky, 1,local_nkx, 1,local_nz)
  CALL allocate_array( inv_ampere_op, 1,local_nky, 1,local_nkx, 1,local_nz)
  CALL allocate_array(   inv_pol_ion, 1,local_nky, 1,local_nkx, 1,local_nz)
  
  ! smaller arrays
  CALL allocate_array(     xnapj, 1,local_na, 1,local_np, 1,local_nj)

  CALL allocate_array(   xnapp2j, 1,local_na, 1,local_np)
  CALL allocate_array(   xnapp1j, 1,local_na, 1,local_np)
  CALL allocate_array(   xnapm1j, 1,local_na, 1,local_np)
  CALL allocate_array(   xnapm2j, 1,local_na, 1,local_np)

  CALL allocate_array(   xnapjp1, 1,local_na, 1,local_nj)
  CALL allocate_array(   xnapjm1, 1,local_na, 1,local_nj)

  CALL allocate_array(   ynapp1j, 1,local_na, 1,local_np, 1,local_nj)
  CALL allocate_array(   ynapm1j, 1,local_na, 1,local_np, 1,local_nj)
  CALL allocate_array( ynapp1jm1, 1,local_na, 1,local_np, 1,local_nj)
  CALL allocate_array( ynapm1jm1, 1,local_na, 1,local_np, 1,local_nj)

  CALL allocate_array(   znapm1j, 1,local_na, 1,local_np, 1,local_nj)
  CALL allocate_array( znapm1jp1, 1,local_na, 1,local_np, 1,local_nj)
  CALL allocate_array( znapm1jm1, 1,local_na, 1,local_np, 1,local_nj)

  CALL allocate_array( xphij,   1,local_na, 1,local_np, 1,local_nj)
  CALL allocate_array( xphijp1, 1,local_na, 1,local_np, 1,local_nj)
  CALL allocate_array( xphijm1, 1,local_na, 1,local_np, 1,local_nj)
  CALL allocate_array( xpsij,   1,local_na, 1,local_np, 1,local_nj)
  CALL allocate_array( xpsijp1, 1,local_na, 1,local_np, 1,local_nj)
  CALL allocate_array( xpsijm1, 1,local_na, 1,local_np, 1,local_nj)
  
  CALL allocate_array( dnjs, 1,jmax+1, 1,jmax+1, 1,jmax+1)
  CALL allocate_array( dv4_Hp_coeff, -2, pmax)
  
  !___________________ 2x5D ARRAYS __________________________
  !! Collision matrices
  IF (GK_CO) THEN !GK collision matrices (one for each kperp)
    CALL allocate_array(  Cab_F, 1,na, 1,na, 1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), 1,local_nky, 1,local_nkx, 1,local_nz)
    CALL allocate_array(  Cab_T, 1,na, 1,na, 1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), 1,local_nky, 1,local_nkx, 1,local_nz)
    CALL allocate_array(  Caa,   1,na,       1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), 1,local_nky, 1,local_nkx, 1,local_nz)
    CALL allocate_array(nuCself, 1,na,       1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), 1,local_nky, 1,local_nkx, 1,local_nz)
  ELSE !DK collision matrix (same for every k)
      CALL allocate_array(  Cab_F, 1,na, 1,na, 1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), 1,1, 1,1, 1,1)
      CALL allocate_array(  Cab_T, 1,na, 1,na, 1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), 1,1, 1,1, 1,1)
      CALL allocate_array(  Caa,    1,na,      1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), 1,1, 1,1, 1,1)
      CALL allocate_array(nuCself,  1,na,      1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), 1,1, 1,1, 1,1)
ENDIF
END SUBROUTINE memory