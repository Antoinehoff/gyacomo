MODULE fields

  use prec_const
  implicit none

!------------------DENSITY Na00------------------

  ! Natural logarithm of normalized density: \ln \hat{N}
!  real(dp), DIMENSION(:,:), ALLOCATABLE :: theta

  ! Natural logarithm of normalized electrons density: \ln \hat{N}
  real(dp), DIMENSION(:,:,:), ALLOCATABLE :: theta_i ! (kr,kz,updatetlevel)

  ! Natural logarithm of normalized ions      density: \ln \hat{N}
  real(dp), DIMENSION(:,:,:), ALLOCATABLE :: theta_e ! (kr,kz,updatetlevel)

!------------------FLUID VELOCITY Na10------------------

  ! Normalized electron gyrocenter parallel fluid velocity \hat{u}
!  real(dp), DIMENSION(:,:), ALLOCATABLE :: vpar

  ! Normalized electron gyrocenter parallel fluid velocity \hat{u}
  real(dp), DIMENSION(:,:,:), ALLOCATABLE :: v_par_e  ! (kr,kz,updatetlevel)

  ! Normalized ion gyrocenter parallel fluid velocity \hat{u}
  real(dp), DIMENSION(:,:,:), ALLOCATABLE :: v_par_i  ! (kr,kz,updatetlevel)

!------------------PARALLEL TEMPERATURE Na20------------------

  ! Natural logarithm of normalized parallel electron temperature: \ln \hat{T} 
! real(dp), DIMENSION(:,:,:), ALLOCATABLE :: temp

  ! Natural logarithm of normalized parallel      electron temperature: \ln \hat{T} 
  real(dp), DIMENSION(:,:,:), ALLOCATABLE :: temp_par_e  ! (kr,kz,updatetlevel)

  ! Natural logarithm of normalized perpendicular electron temperature: \ln \hat{T} 
  real(dp), DIMENSION(:,:,:), ALLOCATABLE :: temp_par_i  ! (kr,kz,updatetlevel)

!------------------PERPENDICULAR TEMPERATURE Na01------------------

  ! Natural logarithm of normalized parallel      ion      temperature: \ln \hat{T} 
  real(dp), DIMENSION(:,:,:), ALLOCATABLE :: temp_per_e  ! (kr,kz,updatetlevel)

  ! Natural logarithm of normalized perpendicular ion      temperature: \ln \hat{T} 
  real(dp), DIMENSION(:,:,:), ALLOCATABLE :: temp_per_i  ! (kr,kz,updatetlevel)

!------------------MOMENTS Napj------------------

  ! Hermite-Moments: N_e^pj ! dimensions correspond to: moments p (from 3 to pmax), z, updatetlevel.
!  real(dp), DIMENSION(:,:,:), ALLOCATABLE :: moments

  ! Hermite-Moments: N_e^pj ! dimensions correspond to: moments pj (from 3 to pmax, 1 to jmax), kr, kz, updatetlevel.
  real(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: moments_e

  ! Hermite-Moments: N_i^pj ! dimensions correspond to: moments pj (from 3 to pmax, 1 to jmax), kr, kz, updatetlevel.
  real(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: moments_i

!------------------ELECTROSTATIC POTENTIAL------------------

  ! Normalized electric potential: \hat{\phi}
!  real(dp), DIMENSION(:,:), ALLOCATABLE :: phi

  ! Normalized electric potential: \hat{\phi}
  real(dp), DIMENSION(:,:,:), ALLOCATABLE :: phi ! (kr,kz,updatetlevel)

END MODULE fields

