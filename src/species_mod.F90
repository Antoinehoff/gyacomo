MODULE species_mod
  USE :: basic
  IMPLICIT NONE
  PRIVATE

  ! Classe that encapsulate all atributes and methods for one arbitrary species
  TYPE, PUBLIC :: species_class
      REAL(dp), PUBLIC :: q     !charge
      REAL(dp), PUBLIC :: sigma !sqrt masse ratio w.r.t. ion mass
      REAL(dp), PUBLIC :: tau   !temperatrue ratio w.r.t. electron temp.
      INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: parray ! Hermite degrees
      INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: jarray ! Laguerre degrees
      ! Hermite-Moments: N_a^pj ! DIMENSIONs correspond to: p, j, kx, ky, z, updatetlevel.
      COMPLEX(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: moments
      ! Arrays to store the rhs, for time integration (ip,ij,ikx,iky,iz,updatetlevel)
      COMPLEX(dp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: moments_rhs
      ! Non linear term array (ip,ij,ikx,iky,iz)
      COMPLEX(dp), DIMENSION(:,:,:,:,:),   ALLOCATABLE :: Sapj
      ! lin rhs p,j coefficient storage (ip,ij)
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xnapj
      REAL(dp), DIMENSION(:),   ALLOCATABLE :: xnapp1j, xnapm1j,   xnapp2j,   xnapm2j, xnapjp1, xnapjm1
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: ynapp1j, ynapm1j,   ynapp1jm1, ynapm1jm1 ! mirror lin coeff for non adiab mom
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: zNapm1j, zNapm1jp1, zNapm1jm1            ! mirror lin coeff for adiab mom
      ! Kernel function evaluation (ij,ikx,iky,iz)
      REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: kernel
      !! Diagnostics
      ! Gyrocenter density (ikx,iky,iz)
      COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: Na00

      ! particle density (ikx,iky,iz)
      COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: density

      ! particle temperature for electron and ions (ikx,iky,iz)
      COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE :: temperature

    CONTAINS
      ! Initialization procedures
      PROCEDURE, PUBLIC :: init                => species_init
      PROCEDURE, PUBLIC :: setup_arrays        => species_setup_arrays
      PROCEDURE, PUBLIC :: evaluate_kernels    => species_evaluate_kernels
      ! Diagnostics
      PROCEDURE, PUBLIC :: compute_density     => species_compute_density
      PROCEDURE, PUBLIC :: compute_temperature => species_compute_temperature
    END TYPE species_class

  ! Routines that every species may use
  CONTAINS
    SUBROUTINE species_setup_arrays(this)

    END SUBROUTINE

    SUBROUTINE species_compute_density(this)

    END SUBROUTINE

END MODULE
