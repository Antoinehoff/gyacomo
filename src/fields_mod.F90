MODULE fields

  use prec_const
  implicit none
!------------------MOMENTS Napj------------------
  ! Hermite-Moments: N_a^pj ! dimensions correspond to: p, j, kr, kz, updatetlevel.
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: moments_e
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: moments_i

!------------------ELECTROSTATIC POTENTIAL------------------

  ! Normalized electric potential: \hat{\phi} ! (kr,kz)
  COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: phi

END MODULE fields
