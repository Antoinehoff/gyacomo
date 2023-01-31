MODULE array

  use prec_const
  implicit none

  ! Arrays to store the rhs, for time integration (ip,ij,iky,ikx,iz,updatetlevel)
  COMPLEX(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: moments_rhs_e
  COMPLEX(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: moments_rhs_i

  ! Arrays of non-adiabatique moments
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: nadiab_moments_e
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: nadiab_moments_i

  ! Derivatives and interpolated moments
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: ddz_nepj
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: interp_nepj
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: ddzND_nepj
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: ddz_nipj
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: interp_nipj
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: ddzND_nipj

  ! Arrays to store special initial modes (semi linear simulation)
  ! Zonal ones (ky=0)
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: moments_e_ZF
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: moments_i_ZF
  COMPLEX(dp), DIMENSION(:,:),     ALLOCATABLE :: phi_ZF
  ! non-zonal modes (kx=0)
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: moments_e_NZ
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: moments_i_NZ
  COMPLEX(dp), DIMENSION(:,:),     ALLOCATABLE :: phi_NZ

  ! Non linear term array (ip,ij,iky,ikx,iz)
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Sepj ! electron
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Sipj ! ion

  ! To load collision matrix (ip,ij,iky,ikx,iz)
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Ceepj, CeipjT
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: CeipjF
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Ciipj, CiepjT
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: CiepjF

  ! Collision term (ip,ij,iky,ikx,iz)
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
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xpsij_e, xpsijp1_e, xpsijm1_e
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xpsij_i, xpsijp1_i, xpsijm1_i
  ! Kernel function evaluation (ij,iky,ikx,iz,odd/even p)
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: kernel_e
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: kernel_i

  ! Poisson operator (iky,ikx,iz)
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: inv_poisson_op
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: inv_ampere_op
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: inv_pol_ion
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: HF_phi_correction_operator

  ! Gyrocenter density for electron and ions (iky,ikx,iz)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Ne00
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Ni00

  ! Kinetic spectrum sum_kx,ky(|Napj(z)|^2), (ip,ij,iz) (should be real)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Nepjz
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Nipjz

  ! particle density for electron and ions (iky,ikx,iz)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: dens_e
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: dens_i

  ! particle fluid velocity for electron and ions (iky,ikx,iz)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: upar_e
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: upar_i
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: uper_e
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: uper_i

  ! particle temperature for electron and ions (iky,ikx,iz)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Tpar_e
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Tpar_i
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Tper_e
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Tper_i
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: temp_e
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: temp_i

END MODULE array
