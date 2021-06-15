MODULE array

  use prec_const
  implicit none

  ! Arrays to store the rhs, for time integration
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: moments_rhs_e ! (ip,ij,ikr,ikz,updatetlevel)
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: moments_rhs_i ! (ip,ij,ikr,ikz,updatetlevel)

  ! To load collision matrix
  REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: Ceepj, CeipjT
  REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: CeipjF
  REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: Ciipj, CiepjT
  REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: CiepjF

  ! Collision term
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: TColl_e, TColl_i
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: TColl_e_local, TColl_i_local

  ! dnjs coefficient storage (in, ij, is)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: dnjs

  ! Kernel function evaluation
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: kernel_e
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: kernel_i

  ! Non linear term array (ip,ij,ikr,ikz)
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: Sepj ! electron
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: Sipj ! ion

  ! Gyrocenter density for electron and ions (meant for 2D output)
  COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: Ne00
  COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: Ni00

  ! particle density for electron and ions (meant for 2D output)
  COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: dens_e
  COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: dens_i

  ! particle temperature for electron and ions (meant for 2D output)
  COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: temp_e
  COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: temp_i

END MODULE array
