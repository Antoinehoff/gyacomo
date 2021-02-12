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

  ! Moments and moments rhs
  CALL allocate_array(     moments_e, ipsg_e,ipeg_e, ijsg_e,ijeg_e, ikrs,ikre, ikzs,ikze, 1,ntimelevel )
  CALL allocate_array(     moments_i, ipsg_i,ipeg_i, ijsg_i,ijeg_i, ikrs,ikre, ikzs,ikze, 1,ntimelevel )
  CALL allocate_array( moments_rhs_e, ips_e,ipe_e, ijs_e,ije_e, ikrs,ikre, ikzs,ikze, 1,ntimelevel )
  CALL allocate_array( moments_rhs_i, ips_i,ipe_i, ijs_i,ije_i, ikrs,ikre, ikzs,ikze, 1,ntimelevel )

  ! Electrostatic potential
  CALL allocate_array(phi, ikrs,ikre, ikzs,ikze)

  ! Electron kernel evaluation from ij=0 to jmax+1 for truncation
  CALL allocate_array(Kernel_e, ijs_e-1,ije_e+1, ikrs,ikre, ikzs,ikze)
  ! Ion kernel evaluation
  CALL allocate_array(Kernel_i, ijs_i-1,ije_i+1, ikrs,ikre, ikzs,ikze)

  ! Collision matrix
  IF (CO .EQ. -1) THEN
    CALL allocate_array(  Ceepj, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1))
    CALL allocate_array( CeipjT, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1))
    CALL allocate_array( CeipjF, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxi+1)*(jmaxi+1))

    CALL allocate_array(  Ciipj, 1,(pmaxi+1)*(jmaxi+1), 1,(pmaxi+1)*(jmaxi+1))
    CALL allocate_array( CiepjT, 1,(pmaxi+1)*(jmaxi+1), 1,(pmaxi+1)*(jmaxi+1))
    CALL allocate_array( CiepjF, 1,(pmaxi+1)*(jmaxi+1), 1,(pmaxe+1)*(jmaxe+1))
  ENDIF

  ! Non linear terms and dnjs table
  CALL allocate_array( Sepj, ips_e,ipe_e, ijs_e,ije_e, ikrs,ikre, ikzs,ikze )
  CALL allocate_array( Sipj, ips_i,ipe_i, ijs_i,ije_i, ikrs,ikre, ikzs,ikze )
  CALL allocate_array( dnjs, 1,maxj+1, 1,maxj+1, 1,maxj+1)

END SUBROUTINE memory
