MODULE array

  use prec_const
  implicit none

  ! Arrays to store the rhs, for time integration (ip,ij,ikx,iky,iz,updatetlevel)
  COMPLEX(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: moments_rhs_e
  COMPLEX(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: moments_rhs_i

  ! Arrays of non-adiabatique moments
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: nadiab_moments_e
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: nadiab_moments_i

  ! Non linear term array (ip,ij,ikx,iky,iz)
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Sepj ! electron
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: Sipj ! ion

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
  REAL(dp), DIMENSION(:),   ALLOCATABLE :: xnepp1j, xnepm1j,   xnepp2j,   xnepm2j, xnepjp1, xnepjm1
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: ynepp1j, ynepm1j,   ynepp1jm1, ynepm1jm1 ! mirror lin coeff for non adiab mom
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: zNepm1j, zNepm1jp1, zNepm1jm1            ! mirror lin coeff for adiab mom
  REAL(dp), DIMENSION(:),   ALLOCATABLE :: xnipp1j, xnipm1j, xnipp2j, xnipm2j, xnipjp1, xnipjm1
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: ynipp1j, ynipm1j,   ynipp1jm1, ynipm1jm1 ! mirror lin coeff for non adiab mom
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: zNipm1j, zNipm1jp1, zNipm1jm1            ! mirror lin coeff for adiab mom
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xphij, xphijp1, xphijm1
  ! Geoemtrical operators
  ! Curvature
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Ckxky  ! dimensions: kx, ky, z
  ! Jacobian
  REAL(dp), DIMENSION(:), ALLOCATABLE :: Jacobian ! dimensions: z
  ! Metric
  REAL(dp), DIMENSION(:), ALLOCATABLE :: gxx, gxy, gyy, gxz, gyz
  ! derivatives of magnetic field strength
  REAL(dp), DIMENSION(:), allocatable :: gradzB  ! dimensions: z
  REAL(dp), DIMENSION(:), allocatable :: gradxB
  ! Relative magnetic field strength
  REAL(dp), DIMENSION(:), allocatable :: hatB
  ! Relative strength of major radius
  REAL(dp), DIMENSION(:), allocatable :: hatR
  ! Geometrical factors
  REAL(dp), DIMENSION(:), allocatable :: Gamma1
  REAL(dp), DIMENSION(:), allocatable :: Gamma2
  REAL(dp), DIMENSION(:), allocatable :: Gamma3
  ! Some geometrical coefficients
  REAL(dp), DIMENSION(:) , allocatable :: gradz_coeff  ! 1 / [ J_{xyz} \hat{B} ]
  ! Kernel function evaluation (ij,ikx,iky,iz)
  REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: kernel_e
  REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: kernel_i
  ! Poisson operator (ikx,iky,iz)
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: inv_poisson_op

  !! Diagnostics
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
