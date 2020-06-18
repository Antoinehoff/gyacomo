SUBROUTINE memory
  ! Allocate arrays (done dynamically otherwise size is unknown)

  USE array
  USE basic
  USE fields
  USE fourier_grid
  USE time_integration  

  use prec_const
  IMPLICIT NONE

  ! Moments and moments rhs
  CALL allocate_array(     moments_e, ips_e,ipe_e, ijs_e,ije_e, ikrs,ikre, ikzs,ikze, 1,ntimelevel )
  CALL allocate_array(     moments_i, ips_i,ipe_i, ijs_i,ije_i, ikrs,ikre, ikzs,ikze, 1,ntimelevel )
  CALL allocate_array( moments_rhs_e, ips_e,ipe_e, ijs_e,ije_e, ikrs,ikre, ikzs,ikze, 1,ntimelevel )
  CALL allocate_array( moments_rhs_i, ips_i,ipe_i, ijs_i,ije_i, ikrs,ikre, ikzs,ikze, 1,ntimelevel )

  ! Electrostatic potential
  CALL allocate_array(phi, ikrs,ikre, ikzs,ikze)

END SUBROUTINE memory
