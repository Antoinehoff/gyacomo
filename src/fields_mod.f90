MODULE fields

  use prec_const
  implicit none
!------------------MOMENTS Napj------------------
  ! Hermite-Moments: N_a^pj ! dimensions correspond to: species (a), p, j, kx, ky, z, updatetlevel.
  COMPLEX(xp), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: moments
!------------------ELECTROSTATIC POTENTIAL------------------

  ! Normalized electric potential: \hat{\phi} ! (kx,ky,z)
  COMPLEX(xp), DIMENSION(:,:,:), ALLOCATABLE :: phi

!------------------MAGNETIC VECTOR
  ! Parallel component (A_\parallel)
  COMPLEX(xp), DIMENSION(:,:,:), ALLOCATABLE :: psi
  ! delta B parallel (not implemeneted)
  ! COMPLEX(xp), DIMENSION(:,:,:), ALLOCATABLE :: dBpar

END MODULE fields
