MODULE array

  use prec_const
  implicit none

  ! Arrays to store the rhs, for time integration (ip,ij,ikx,iky,iz,updatetlevel)
  COMPLEX(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: moments_rhs_e
  COMPLEX(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: moments_rhs_i

  ! To load collision matrix (ip1,ij1,ip2,ij2)
  REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: Ceepj, CeipjT
  REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: CeipjF
  REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: Ciipj, CiepjT
  REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: CiepjF

  ! Collision term (ip,ij,ikx,iky,iz)
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: TColl_e, TColl_i
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: TColl_e_local, TColl_i_local

  ! dnjs coefficient storage (in, ij, is)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: dnjs

  ! lin rhs p,j coefficient storage (ip,ij)
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xnepj,xnipj
  REAL(dp), DIMENSION(:),   ALLOCATABLE :: xnepp1j, xnepm1j, xnepp2j, xnepm2j, xnepjp1, xnepjm1
  REAL(dp), DIMENSION(:),   ALLOCATABLE :: xnipp1j, xnipm1j, xnipp2j, xnipm2j, xnipjp1, xnipjm1
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xphij, xphijp1, xphijm1

  ! Curvature array
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Ckxky

  ! Kernel function evaluation (ij,ikx,iky)
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: kernel_e
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: kernel_i

  ! Non linear term array (ip,ij,ikx,iky,iz)
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Sepj ! electron
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Sipj ! ion

  ! Gyrocenter density for electron and ions (ikx,iky,iz)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Ne00
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Ni00

  ! particle density for electron and ions (ikx,iky,iz)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: dens_e
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: dens_i

  ! particle temperature for electron and ions (ikx,iky,iz)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: temp_e
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: temp_i

END MODULE array
