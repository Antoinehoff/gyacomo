MODULE array

  use prec_const
  implicit none

  ! Arrays to store the rhs, for time integration (ia,ip,ij,iky,ikx,iz)
  COMPLEX(dp), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: moments_rhs

  ! Arrays of non-adiabatique moments
  COMPLEX(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: nadiab_moments

  ! Derivatives and interpolated moments
  COMPLEX(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: ddz_napj
  COMPLEX(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: interp_napj
  COMPLEX(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: ddzND_Napj

  ! Non linear term array (ip,ij,iky,ikx,iz)
  COMPLEX(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: Sapj ! electron

  ! a-a collision matrix (ia,ip,ij,iky,ikx,iz)
  REAL(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: Caa
  ! Test and field collision matrices (ia,ib,ip,ij,iky,ikx,iz)
  REAL(dp), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: Cab_F, Cab_T
  ! nu x self collision matrix nuCself = nuaa*Caa + sum_b_neq_a nu_ab Cab_T
  REAL(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: nuCself

  ! Collision term (ip,ij,iky,ikx,iz)
  COMPLEX(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: Capj
  COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE   :: TColl_e_local, TColl_i_local

  ! dnjs coefficient storage (in, ij, is)
  COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: dnjs

  ! Hermite fourth derivative coeff storage 4*sqrt(p!/(p-4)!)
  COMPLEX(dp), DIMENSION(:), ALLOCATABLE :: dv4_Hp_coeff

  ! lin rhs p,j coefficient storage (ip,ij)
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: xnapj
  REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: xnapp1j, xnapp2j,   xnapm1j, xnapm2j, xnapjp1, xnapjm1
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ynapp1j, ynapm1j,   ynapp1jm1, ynapm1jm1 ! mirror lin coeff for non adiab mom
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: zNapm1j, zNapm1jp1, zNapm1jm1            ! mirror lin coeff for adiab mom
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: xphij, xphijp1, xphijm1
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: xpsij, xpsijp1, xpsijm1
  ! Kernel function evaluation (ia,ij,iky,ikx,iz,odd/even p)
  REAL(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: kernel

  ! Poisson operator (iky,ikx,iz)
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: inv_poisson_op
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: inv_ampere_op
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: inv_pol_ion
  ! REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: HF_phi_correction_operator

  ! Kinetic spectrum sum_kx,ky(|Napj(z)|^2), (ia,ip,ij,iz) (should be real)
  REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: Napjz

  ! particle density for electron and ions (iky,ikx,iz)
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: dens

  ! particle fluid velocity for electron and ions (iky,ikx,iz)
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: upar
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: uper

  ! particle temperature for electron and ions (iky,ikx,iz)
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: Tpar
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: Tper
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: temp

END MODULE array
