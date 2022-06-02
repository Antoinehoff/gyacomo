SUBROUTINE memory
  ! Allocate arrays (done dynamically otherwise size is unknown)

  USE array
  USE basic
  USE fields
  USE grid
  USE time_integration
  USE model, ONLY: LINEARITY, KIN_E
  USE collision

  USE prec_const
  IMPLICIT NONE

  ! Electrostatic potential
  CALL allocate_array(           phi, ikys,ikye, ikxs,ikxe, izgs,izge)
  CALL allocate_array(        phi_ZF, ikxs,ikxe, izs,ize)
  CALL allocate_array(        phi_EM, ikys,ikye, izs,ize)
  CALL allocate_array(inv_poisson_op, ikys,ikye, ikxs,ikxe, izs,ize)

  !Electrons arrays
  IF(KIN_E) THEN
  CALL allocate_array(             Ne00, ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(           dens_e, ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(           upar_e, ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(           uper_e, ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(           Tpar_e, ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(           Tper_e, ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(           temp_e, ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(         Kernel_e,                ijgs_e,ijge_e, ikys,ikye, ikxs,ikxe, izgs,izge,  0,1)
  CALL allocate_array(        moments_e, ipgs_e,ipge_e, ijgs_e,ijge_e, ikys,ikye, ikxs,ikxe, izgs,izge, 1,ntimelevel )
  CALL allocate_array(            Nepjz,  ips_e,ipe_e,   ijs_e,ije_e,                         izs,ize)
  CALL allocate_array(    moments_rhs_e,  ips_e,ipe_e,   ijs_e,ije_e,  ikys,ikye, ikxs,ikxe,  izs,ize,  1,ntimelevel )
  CALL allocate_array( nadiab_moments_e, ipgs_e,ipge_e, ijgs_e,ijge_e, ikys,ikye, ikxs,ikxe, izgs,izge)
  CALL allocate_array(         ddz_nepj, ipgs_e,ipge_e, ijgs_e,ijge_e, ikys,ikye, ikxs,ikxe, izgs,izge)
  CALL allocate_array(        ddz4_Nepj, ipgs_e,ipge_e, ijgs_e,ijge_e, ikys,ikye, ikxs,ikxe, izgs,izge)
  CALL allocate_array(      interp_nepj, ipgs_e,ipge_e, ijgs_e,ijge_e, ikys,ikye, ikxs,ikxe, izgs,izge)
  CALL allocate_array(     moments_e_ZF, ipgs_e,ipge_e, ijgs_e,ijge_e, ikxs,ikxe, izs,ize)
  CALL allocate_array(     moments_e_EM, ipgs_e,ipge_e, ijgs_e,ijge_e, ikys,ikye, izs,ize)
  CALL allocate_array(          TColl_e,  ips_e,ipe_e,   ijs_e,ije_e , ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(             Sepj,  ips_e,ipe_e,   ijs_e,ije_e,  ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(           xnepj,   ips_e,ipe_e,   ijs_e,ije_e)
  CALL allocate_array(           xnepp2j, ips_e,ipe_e)
  CALL allocate_array(           xnepp1j, ips_e,ipe_e)
  CALL allocate_array(           xnepm1j, ips_e,ipe_e)
  CALL allocate_array(           xnepm2j, ips_e,ipe_e)
  CALL allocate_array(           xnepjp1,                ijs_e,ije_e)
  CALL allocate_array(           xnepjm1,                ijs_e,ije_e)
  CALL allocate_array(           ynepp1j, ips_e,ipe_e,   ijs_e,ije_e)
  CALL allocate_array(           ynepm1j, ips_e,ipe_e,   ijs_e,ije_e)
  CALL allocate_array(         ynepp1jm1, ips_e,ipe_e,   ijs_e,ije_e)
  CALL allocate_array(         ynepm1jm1, ips_e,ipe_e,   ijs_e,ije_e)
  CALL allocate_array(           zNepm1j, ips_e,ipe_e,   ijs_e,ije_e)
  CALL allocate_array(         zNepm1jp1, ips_e,ipe_e,   ijs_e,ije_e)
  CALL allocate_array(         zNepm1jm1, ips_e,ipe_e,   ijs_e,ije_e)
  ENDIF

  !Ions arrays
  CALL allocate_array(             Ni00, ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(           dens_i, ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(           upar_i, ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(           uper_i, ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(           Tpar_i, ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(           Tper_i, ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(           temp_i, ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(         Kernel_i,                ijgs_i,ijge_i, ikys,ikye, ikxs,ikxe, izgs,izge,  0,1)
  CALL allocate_array(        moments_i, ipgs_i,ipge_i, ijgs_i,ijge_i, ikys,ikye, ikxs,ikxe, izgs,izge, 1,ntimelevel )
  CALL allocate_array(            Nipjz,  ips_i,ipe_i,   ijs_i,ije_i,                         izs,ize)
  CALL allocate_array(    moments_rhs_i,  ips_i,ipe_i,   ijs_i,ije_i,  ikys,ikye, ikxs,ikxe,  izs,ize,  1,ntimelevel )
  CALL allocate_array( nadiab_moments_i, ipgs_i,ipge_i, ijgs_i,ijge_i, ikys,ikye, ikxs,ikxe, izgs,izge)
  CALL allocate_array(         ddz_nipj, ipgs_i,ipge_i, ijgs_i,ijge_i, ikys,ikye, ikxs,ikxe, izgs,izge)
  CALL allocate_array(        ddz4_Nipj, ipgs_i,ipge_i, ijgs_i,ijge_i, ikys,ikye, ikxs,ikxe, izgs,izge)
  CALL allocate_array(      interp_nipj, ipgs_i,ipge_i, ijgs_i,ijge_i, ikys,ikye, ikxs,ikxe, izgs,izge)
  CALL allocate_array(     moments_i_ZF, ipgs_i,ipge_i, ijgs_i,ijge_i, ikxs,ikxe, izs,ize)
  CALL allocate_array(     moments_i_EM, ipgs_i,ipge_i, ijgs_i,ijge_i, ikys,ikye, izs,ize)
  CALL allocate_array(          TColl_i,  ips_i,ipe_i,   ijs_i,ije_i,  ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array(             Sipj,  ips_i,ipe_i,   ijs_i,ije_i,  ikys,ikye, ikxs,ikxe, izs,ize)
  CALL allocate_array( xnipj,   ips_i,ipe_i, ijs_i,ije_i)
  CALL allocate_array( xnipp2j, ips_i,ipe_i)
  CALL allocate_array( xnipp1j, ips_i,ipe_i)
  CALL allocate_array( xnipm1j, ips_i,ipe_i)
  CALL allocate_array( xnipm2j, ips_i,ipe_i)
  CALL allocate_array( xnipjp1, ijs_i,ije_i)
  CALL allocate_array( xnipjm1, ijs_i,ije_i)
  CALL allocate_array(   ynipp1j, ips_i,ipe_i, ijs_i,ije_i)
  CALL allocate_array(   ynipm1j, ips_i,ipe_i, ijs_i,ije_i)
  CALL allocate_array( ynipp1jm1, ips_i,ipe_i, ijs_i,ije_i)
  CALL allocate_array( ynipm1jm1, ips_i,ipe_i, ijs_i,ije_i)
  CALL allocate_array(   zNipm1j, ips_i,ipe_i, ijs_i,ije_i)
  CALL allocate_array( zNipm1jp1, ips_i,ipe_i, ijs_i,ije_i)
  CALL allocate_array( zNipm1jm1, ips_i,ipe_i, ijs_i,ije_i)

  ! Non linear terms and dnjs table
  CALL allocate_array( dnjs, 1,maxj+1, 1,maxj+1, 1,maxj+1)

  ! elect. pot. linear terms
  IF (KIN_E) THEN
    CALL allocate_array( xphij_e,   ips_e,ipe_e, ijs_e,ije_e)
    CALL allocate_array( xphijp1_e, ips_e,ipe_e, ijs_e,ije_e)
    CALL allocate_array( xphijm1_e, ips_e,ipe_e, ijs_e,ije_e)
  ENDIF
  CALL allocate_array( xphij_i,   ips_i,ipe_i, ijs_i,ije_i)
  CALL allocate_array( xphijp1_i, ips_i,ipe_i, ijs_i,ije_i)
  CALL allocate_array( xphijm1_i, ips_i,ipe_i, ijs_i,ije_i)

  !___________________ 2x5D ARRAYS __________________________
  !! Collision matrices
  IF (gyrokin_CO) THEN !GK collision matrices (one for each kperp)
    IF (KIN_E) THEN
    CALL allocate_array(  Ceepj, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), ikys,ikye, ikxs,ikxe, izs,ize)
    CALL allocate_array( CeipjT, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), ikys,ikye, ikxs,ikxe, izs,ize)
    CALL allocate_array( CeipjF, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxi+1)*(jmaxi+1), ikys,ikye, ikxs,ikxe, izs,ize)
    CALL allocate_array( CiepjT, 1,(pmaxi+1)*(jmaxi+1), 1,(pmaxi+1)*(jmaxi+1), ikys,ikye, ikxs,ikxe, izs,ize)
    CALL allocate_array( CiepjF, 1,(pmaxi+1)*(jmaxi+1), 1,(pmaxe+1)*(jmaxe+1), ikys,ikye, ikxs,ikxe, izs,ize)
    ENDIF
    CALL allocate_array(  Ciipj, 1,(pmaxi+1)*(jmaxi+1), 1,(pmaxi+1)*(jmaxi+1), ikys,ikye, ikxs,ikxe, izs,ize)
  ELSE !DK collision matrix (same for every k)
      IF (KIN_E) THEN
      CALL allocate_array(  Ceepj, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), 1,1, 1,1, 1,1)
      CALL allocate_array( CeipjT, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), 1,1, 1,1, 1,1)
      CALL allocate_array( CeipjF, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxi+1)*(jmaxi+1), 1,1, 1,1, 1,1)
      CALL allocate_array( CiepjT, 1,(pmaxi+1)*(jmaxi+1), 1,(pmaxi+1)*(jmaxi+1), 1,1, 1,1, 1,1)
      CALL allocate_array( CiepjF, 1,(pmaxi+1)*(jmaxi+1), 1,(pmaxe+1)*(jmaxe+1), 1,1, 1,1, 1,1)
      ENDIF
      CALL allocate_array(  Ciipj, 1,(pmaxi+1)*(jmaxi+1), 1,(pmaxi+1)*(jmaxi+1), 1,1, 1,1, 1,1)
 ENDIF

END SUBROUTINE memory
