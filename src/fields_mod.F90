MODULE fields

  use prec_const
  implicit none
!------------------MOMENTS Napj------------------
  ! Hermite-Moments: N_a^pj ! dimensions correspond to: ipj, kr, kz, updatetlevel.
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: moments

!------------------ELECTROSTATIC POTENTIAL------------------

  ! Normalized electric potential: \hat{\phi} ! (kr,kz)
  COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: phi 

END MODULE fields

