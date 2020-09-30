MODULE array

  use prec_const
  implicit none

  ! Arrays to store the rhs, for time integration
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: moments_rhs_e ! (ip,ij,ikr,ikz,updatetlevel)
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: moments_rhs_i ! (ip,ij,ikr,ikz,updatetlevel)

  ! To load collision matrix
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Ceepj, CeipjT
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: CeipjF
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Ciipj, CiepjT
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: CiepjF

  ! dnjs coefficient storage (in, ij, is)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: dnjs

  ! Non linear term array (ip,ij,ikr,ikz)
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: Sepj ! electron
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: Sipj ! ion

END MODULE array
