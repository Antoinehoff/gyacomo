SUBROUTINE memory
  ! Allocate arrays (done dynamically otherwise size is unknown)

  USE array
  USE basic
  USE fields
  USE grid
  USE time_integration
  USE model, ONLY: CO, NON_LIN

  USE prec_const
  IMPLICIT NONE
  !___________________ 2D ARRAYS __________________________
  ! Electrostatic potential
  CALL allocate_array(phi, ikxs,ikxe, ikys,ikye, izs,ize)

  !! Diagnostics arrays
  ! Gyrocenter density *for 2D output*
  CALL allocate_array(Ne00, ikxs,ikxe, ikys,ikye, izs,ize)
  CALL allocate_array(Ni00, ikxs,ikxe, ikys,ikye, izs,ize)
  ! particle density *for 2D output*
  CALL allocate_array(dens_e, ikxs,ikxe, ikys,ikye, izs,ize)
  CALL allocate_array(dens_i, ikxs,ikxe, ikys,ikye, izs,ize)
  ! particle temperature *for 2D output*
  CALL allocate_array(temp_e, ikxs,ikxe, ikys,ikye, izs,ize)
  CALL allocate_array(temp_i, ikxs,ikxe, ikys,ikye, izs,ize)

  !___________________ 3D ARRAYS __________________________
  !! FLR kernels functions
  ! Kernel evaluation from j= -1 to jmax+1 for truncation
  CALL allocate_array(Kernel_e, ijsg_e,ijeg_e, ikxs,ikxe, ikys,ikye)
  CALL allocate_array(Kernel_i, ijsg_i,ijeg_i, ikxs,ikxe, ikys,ikye)

  !___________________ 5D ARRAYS __________________________
  ! Moments with ghost degrees for p+2 p-2, j+1, j-1 truncations
  CALL allocate_array( moments_e, ipsg_e,ipeg_e, ijsg_e,ijeg_e, ikxs,ikxe, ikys,ikye, izs,ize, 1,ntimelevel )
  CALL allocate_array( moments_i, ipsg_i,ipeg_i, ijsg_i,ijeg_i, ikxs,ikxe, ikys,ikye, izs,ize, 1,ntimelevel )

  ! Moments right-hand-side (contains linear part of hierarchy)
  CALL allocate_array( moments_rhs_e, ips_e,ipe_e, ijs_e,ije_e, ikxs,ikxe, ikys,ikye, izs,ize, 1,ntimelevel )
  CALL allocate_array( moments_rhs_i, ips_i,ipe_i, ijs_i,ije_i, ikxs,ikxe, ikys,ikye, izs,ize, 1,ntimelevel )

  ! Collision term
  CALL allocate_array(  TColl_e, ips_e,ipe_e, ijs_e,ije_e , ikxs,ikxe, ikys,ikye, izs,ize)
  CALL allocate_array(  TColl_i, ips_i,ipe_i, ijs_i,ije_i , ikxs,ikxe, ikys,ikye, izs,ize)

  ! Non linear terms and dnjs table
  CALL allocate_array( Sepj, ips_e,ipe_e, ijs_e,ije_e, ikxs,ikxe, ikys,ikye, izs,ize)
  CALL allocate_array( Sipj, ips_i,ipe_i, ijs_i,ije_i, ikxs,ikxe, ikys,ikye, izs,ize)
  CALL allocate_array( dnjs, 1,maxj+1, 1,maxj+1, 1,maxj+1)

  ! Linear coeff for moments rhs
  ! electrons
  CALL allocate_array( xnepj,   ips_e,ipe_e, ijs_e,ije_e)
  CALL allocate_array( xnepp2j, ips_e,ipe_e)
  CALL allocate_array( xnepp1j, ips_e,ipe_e)
  CALL allocate_array( xnepm1j, ips_e,ipe_e)
  CALL allocate_array( xnepm2j, ips_e,ipe_e)
  CALL allocate_array( xnepjp1, ijs_e,ije_e)
  CALL allocate_array( xnepjm1, ijs_e,ije_e)
  ! ions
  CALL allocate_array( xnipj,   ips_i,ipe_i, ijs_i,ije_i)
  CALL allocate_array( xnipp2j, ips_i,ipe_i)
  CALL allocate_array( xnipp1j, ips_i,ipe_i)
  CALL allocate_array( xnipm1j, ips_i,ipe_i)
  CALL allocate_array( xnipm2j, ips_i,ipe_i)
  CALL allocate_array( xnipjp1, ijs_i,ije_i)
  CALL allocate_array( xnipjm1, ijs_i,ije_i)
  ! elect. pot.
  CALL allocate_array( xphij,   ips_i,ipe_i, ijs_i,ije_i)
  CALL allocate_array( xphijp1, ips_i,ipe_i, ijs_i,ije_i)
  CALL allocate_array( xphijm1, ips_i,ipe_i, ijs_i,ije_i)

  ! Curvature and geometry
  CALL allocate_array( Ckxky, ikxs,ikxe, ikys,ikye, izs,ize)

  !___________________ 2x5D ARRAYS __________________________
  !! Collision matrices
  IF (CO .LT. -1) THEN !DK collision matrix (same for every k)
    CALL allocate_array(  Ceepj, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), 1,1, 1,1)
    CALL allocate_array( CeipjT, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), 1,1, 1,1)
    CALL allocate_array( CeipjF, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxi+1)*(jmaxi+1), 1,1, 1,1)
    CALL allocate_array(  Ciipj, 1,(pmaxi+1)*(jmaxi+1), 1,(pmaxi+1)*(jmaxi+1), 1,1, 1,1)
    CALL allocate_array( CiepjT, 1,(pmaxi+1)*(jmaxi+1), 1,(pmaxi+1)*(jmaxi+1), 1,1, 1,1)
    CALL allocate_array( CiepjF, 1,(pmaxi+1)*(jmaxi+1), 1,(pmaxe+1)*(jmaxe+1), 1,1, 1,1)
  ELSEIF (CO .GT. 1) THEN !GK collision matrices (one for each kperp)
    CALL allocate_array(  Ceepj, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), ikxs,ikxe, ikys,ikye)
    CALL allocate_array( CeipjT, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), ikxs,ikxe, ikys,ikye)
    CALL allocate_array( CeipjF, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxi+1)*(jmaxi+1), ikxs,ikxe, ikys,ikye)
    CALL allocate_array(  Ciipj, 1,(pmaxi+1)*(jmaxi+1), 1,(pmaxi+1)*(jmaxi+1), ikxs,ikxe, ikys,ikye)
    CALL allocate_array( CiepjT, 1,(pmaxi+1)*(jmaxi+1), 1,(pmaxi+1)*(jmaxi+1), ikxs,ikxe, ikys,ikye)
    CALL allocate_array( CiepjF, 1,(pmaxi+1)*(jmaxi+1), 1,(pmaxe+1)*(jmaxe+1), ikxs,ikxe, ikys,ikye)
  ENDIF

END SUBROUTINE memory
