MODULE fields

  use prec_const
  implicit none
!------------------MOMENTS Napj------------------
  ! Hermite-Moments: N_a^pj ! dimensions correspond to: species (a), p, j, kx, ky, z, updatetlevel.
  COMPLEX(xp), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: moments
!------------------ELECTROSTATIC POTENTIAL------------------

  ! Normalized electric potential: \hat{\phi} ! (kx,ky,z)
  COMPLEX(xp), DIMENSION(:,:,:), ALLOCATABLE :: phi

!------------------Vector field part
  COMPLEX(xp), DIMENSION(:,:,:), ALLOCATABLE :: psi

END MODULE fields
