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
  ! Parallel magnetic fluctuation, normalized to the local field:
  !   dBpar = (delta B_par/B0(z)) * (L_ref/rho_s)          ! (kx,ky,z)
  ! Solved with phi from perpendicular pressure balance when BPAR is on,
  ! identically zero otherwise.
  COMPLEX(xp), DIMENSION(:,:,:), ALLOCATABLE :: dBpar

END MODULE fields
