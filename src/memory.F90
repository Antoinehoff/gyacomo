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
  CALL allocate_array(     moments, ipjs,ipje, ikrs,ikre, ikzs,ikze, 1,ntimelevel)
  CALL allocate_array( moments_rhs, ipjs,ipje, ikrs,ikre, ikzs,ikze, 1,ntimelevel)

  ! Electrostatic potential
  CALL allocate_array(         phi, ikrs,ikre, ikzs,ikze)

  ! Intermediate steps in rhs of equations
  !CALL allocate_array( moments_Apl, ipjs,ipje, ikrs,ikre, ikzs,ikze,)
  !CALL allocate_array( moments_Bpl, ipjs,ipje, ikrs,ikre, ikzs,ikze,)
  !CALL allocate_array( moments_Cpl, ipjs,ipje, ikrs,ikre, ikzs,ikze,)
  !CALL allocate_array( moments_Dpl, ipjs,ipje, ikrs,ikre, ikzs,ikze,)
  !CALL allocate_array( moments_Epl, ipjs,ipje, ikrs,ikre, ikzs,ikze,)
  !CALL allocate_array( moments_Fpl, ipjs,ipje, ikrs,ikre, ikzs,ikze,)
  !CALL allocate_array( moments_Gpl, ipjs,ipje, ikrs,ikre, ikzs,ikze,)
  !CALL allocate_array( moments_Hpl, ipjs,ipje, ikrs,ikre, ikzs,ikze,)
END SUBROUTINE memory
