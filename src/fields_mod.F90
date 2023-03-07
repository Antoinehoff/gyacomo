MODULE fields

  use prec_const
  implicit none
!------------------MOMENTS Napj------------------
  ! Hermite-Moments: N_a^pj ! dimensions correspond to: species (a), p, j, kx, ky, z, updatetlevel.
  COMPLEX(dp), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: moments
!------------------ELECTROSTATIC POTENTIAL------------------

  ! Normalized electric potential: \hat{\phi} ! (kx,ky,z)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: phi

!------------------Vector field part
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: psi

END MODULE fields
