MODULE array

  use prec_const
  implicit none

  ! Arrays to store the rhs, for time integration (ia,ip,ij,iky,ikx,iz)
  COMPLEX(xp), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: moments_rhs

  ! Arrays of non-adiabatique moments
  COMPLEX(xp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: nadiab_moments

  ! Derivatives and interpolated moments
  COMPLEX(xp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: ddz_napj
  COMPLEX(xp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: interp_napj
  COMPLEX(xp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: ddzND_Napj

  ! Non linear term array (ip,ij,iky,ikx,iz)
  COMPLEX(xp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: Sapj ! electron

  ! a-a collision matrix (ia,ip,ij,iky,ikx,iz)
  REAL(xp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: Caa
  ! Test and field collision matrices (ia,ib,ip,ij,iky,ikx,iz)
  REAL(xp), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: Cab_F, Cab_T
  ! nu x self collision matrix nuCself = nuaa*Caa + sum_b_neq_a nu_ab Cab_T
  REAL(xp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: nuCself

  ! Collision term (ip,ij,iky,ikx,iz)
  COMPLEX(xp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: Capj
  COMPLEX(xp), DIMENSION(:,:,:,:,:), ALLOCATABLE   :: TColl_e_local, TColl_i_local

  ! dnjs coefficient storage (in, ij, is)
  COMPLEX(xp), DIMENSION(:,:,:), ALLOCATABLE :: dnjs

  ! Hermite fourth derivative coeff storage 4*sqrt(p!/(p-4)!)
  COMPLEX(xp), DIMENSION(:), ALLOCATABLE :: dv4_Hp_coeff

  ! lin rhs p,j coefficient storage (ip,ij)
  REAL(xp), DIMENSION(:,:,:), ALLOCATABLE :: xnapj
  REAL(xp), DIMENSION(:,:),   ALLOCATABLE :: xnapp1j, xnapp2j,   xnapm1j, xnapm2j, xnapjp1, xnapjm1
  REAL(xp), DIMENSION(:,:,:), ALLOCATABLE :: ynapp1j, ynapm1j,   ynapp1jm1, ynapm1jm1 ! mirror lin coeff for non adiab mom
  REAL(xp), DIMENSION(:,:,:), ALLOCATABLE :: zNapm1j, zNapm1jp1, zNapm1jm1            ! mirror lin coeff for adiab mom
  REAL(xp), DIMENSION(:,:,:), ALLOCATABLE :: xphij, xphijp1, xphijm1
  REAL(xp), DIMENSION(:,:,:), ALLOCATABLE :: xpsij, xpsijp1, xpsijm1
  ! Kernel function evaluation (ia,ij,iky,ikx,iz,odd/even p)
  REAL(xp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: kernel

  ! Poisson operator (iky,ikx,iz)
  REAL(xp), DIMENSION(:,:,:), ALLOCATABLE :: inv_poisson_op
  REAL(xp), DIMENSION(:,:,:), ALLOCATABLE :: inv_ampere_op
  REAL(xp), DIMENSION(:,:,:), ALLOCATABLE :: inv_pol_ion
  ! REAL(xp), DIMENSION(:,:,:), ALLOCATABLE :: HF_phi_correction_operator

  ! Kinetic spectrum sum_kx,ky(|Napj(z)|^2), (ia,ip,ij,iz) (should be real)
  REAL(xp), DIMENSION(:,:,:,:), ALLOCATABLE :: Napjz

  ! particle density for electron and ions (iky,ikx,iz)
  COMPLEX(xp), DIMENSION(:,:,:,:), ALLOCATABLE :: dens

  ! particle fluid velocity for electron and ions (iky,ikx,iz)
  COMPLEX(xp), DIMENSION(:,:,:,:), ALLOCATABLE :: upar
  COMPLEX(xp), DIMENSION(:,:,:,:), ALLOCATABLE :: uper

  ! particle temperature for electron and ions (iky,ikx,iz)
  COMPLEX(xp), DIMENSION(:,:,:,:), ALLOCATABLE :: Tpar
  COMPLEX(xp), DIMENSION(:,:,:,:), ALLOCATABLE :: Tper
  COMPLEX(xp), DIMENSION(:,:,:,:), ALLOCATABLE :: temp

END MODULE array
