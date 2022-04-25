MODULE array

  use prec_const
  implicit none

  ! Arrays to store the rhs, for time integration (ip,ij,ikx,iky,iz,updatetlevel)
  COMPLEX(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: moments_rhs_e
  COMPLEX(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: moments_rhs_i

  ! Arrays of non-adiabatique moments
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: nadiab_moments_e
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: nadiab_moments_i

  ! Arrays to store special initial modes (semi linear simulation)
  ! Zonal ones (ky=0)
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: moments_e_ZF
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: moments_i_ZF
  COMPLEX(dp), DIMENSION(:,:),     ALLOCATABLE :: phi_ZF
  ! Entropy modes (kx=0)
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: moments_e_EM
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: moments_i_EM
  COMPLEX(dp), DIMENSION(:,:),     ALLOCATABLE :: phi_EM

  ! Non linear term array (ip,ij,ikx,iky,iz)
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Sepj ! electron
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Sipj ! ion

  ! To load collision matrix (ip,ij,ikx,iky,iz)
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Ceepj, CeipjT
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: CeipjF
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Ciipj, CiepjT
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: CiepjF

  ! Collision term (ip,ij,ikx,iky,iz)
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: TColl_e, TColl_i
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: TColl_e_local, TColl_i_local

  ! dnjs coefficient storage (in, ij, is)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: dnjs

  ! lin rhs p,j coefficient storage (ip,ij)
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xnepj,xnipj
  REAL(dp), DIMENSION(:),   ALLOCATABLE :: xnepp1j, xnepm1j,   xnepp2j,   xnepm2j, xnepjp1, xnepjm1
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: ynepp1j, ynepm1j,   ynepp1jm1, ynepm1jm1 ! mirror lin coeff for non adiab mom
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: zNepm1j, zNepm1jp1, zNepm1jm1            ! mirror lin coeff for adiab mom
  REAL(dp), DIMENSION(:),   ALLOCATABLE :: xnipp1j, xnipm1j, xnipp2j, xnipm2j, xnipjp1, xnipjm1
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: ynipp1j, ynipm1j,   ynipp1jm1, ynipm1jm1 ! mirror lin coeff for non adiab mom
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: zNipm1j, zNipm1jp1, zNipm1jm1            ! mirror lin coeff for adiab mom
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xphij_e, xphijp1_e, xphijm1_e
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xphij_i, xphijp1_i, xphijm1_i
  
  ! Kernel function evaluation (ij,ikx,iky,iz,odd/even p)
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: kernel_e
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: kernel_i
  ! Poisson operator (ikx,iky,iz)
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: inv_poisson_op

  !! Diagnostics
  ! Gyrocenter density for electron and ions (ikx,iky,iz)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Ne00
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Ni00

  ! particle density for electron and ions (ikx,iky,iz)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: dens_e
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: dens_i

  ! particle fluid velocity for electron and ions (ikx,iky,iz)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: upar_e
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: upar_i
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: uper_e
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: uper_i

  ! particle temperature for electron and ions (ikx,iky,iz)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Tpar_e
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Tpar_i
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Tper_e
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Tper_i
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: temp_e
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: temp_i

END MODULE array
