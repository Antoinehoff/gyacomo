MODULE array

  !USE mumps_bsplines
  use prec_const
  implicit none

  ! Fields, on opposite grid
  real(dp), DIMENSION(:), ALLOCATABLE :: vpar_n

  ! Arrays to store the rhs, for time integration
  real(dp), DIMENSION(:,:,:), ALLOCATABLE :: theta_rhs, temp_rhs, vpar_rhs ! (ikr,ikz,updatetlevel)
  real(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: moments_rhs ! (ip,ij,ikr,ikz,updatetlevel)
  ! rhs, on opposite grid
  real(dp), DIMENSION(:,:), ALLOCATABLE :: temp_rhs_v, vpar_rhs_n !(ikr,ikz)


  real(dp), DIMENSION(:,:), ALLOCATABLE:: sqrt_exp_temp ! Often needed in 'rhs', compute it only once and store it
  real(dp), DIMENSION(:,:), ALLOCATABLE:: sqrt_exp_temp_v ! On opposite grid

  ! Spatial 1rst derivatives on respective grids
  ! real(dp), DIMENSION(:,:), ALLOCATABLE:: thetaz, sqrt_exp_tempz, vparz, phiz ! (ikr,ikz)
  ! real(dp), DIMENSION(:,:,:,:), ALLOCATABLE:: momentsz ! (ip,ij,ikr,ikz)
  ! Spatial 1rst derivatives on opposite grid
  ! real(dp), DIMENSION(:,:), ALLOCATABLE:: thetaz_v, sqrt_exp_tempz_v, vparz_n, phiz_v ! (ikr,ikz)
  
  ! ! For upwind scheme, other set of spatial derivatives, the "minus side"
  ! real(dp), DIMENSION(:), ALLOCATABLE:: thetazm, tempzm, vparzm, phizm
  ! real(dp), DIMENSION(:,:), ALLOCATABLE:: momentszm

  ! Spatial 2nd derivatives, for numerical diffusion, on respective grids
  ! real(dp), DIMENSION(:,:), ALLOCATABLE:: thetazz, tempzz, vparzz
  ! real(dp), DIMENSION(:,:,:), ALLOCATABLE:: momentszz

  ! Intermediate steps in rhs of equations
  real(dp), DIMENSION(:,:), ALLOCATABLE:: theta_Cpl, theta_Hpl
  real(dp), DIMENSION(:,:), ALLOCATABLE:: temp_Cpl, temp_Fpl, temp_Hpl, temp_Ipl
  real(dp), DIMENSION(:,:), ALLOCATABLE:: vpar_Cpl, vpar_Dpl, vpar_Epl, vpar_Fpl, vpar_Hpl
  real(dp), DIMENSION(:,:), ALLOCATABLE:: moments_Apl, moments_Bpl, moments_Cpl, moments_Dpl, moments_Epl, moments_Fpl, moments_Gpl, moments_Hpl
  ! moments_Apl = A^l
  ! moments_Bpl = sum_p [ B_p^l N_e^p ]
  ! moments_Cpl = sum_p [ C_p^l N_e^p ]
  ! moments_Dpl = sum_p [ D_p^l N_e^p ]
  ! moments_Epl = sum_p [ E_p^l N_e^p ]
  ! moments_Fpl = sum_p [ F_p^l N_e^p ]
  ! moments_Gpl = sum_p [ G_p^l N_e^p ]
  ! moments_Hpl = sum_p [ H_p^l N_e^p ]
  real(dp) :: moments_Ipl
  ! moments_Ipl = sum_p [ I_p^l d(N_e^p)/dz ]


  ! For poisson solver
  !TYPE(mumps_mat), SAVE :: laplace_smumps
  real(dp), DIMENSION(:,:), ALLOCATABLE :: laplace_rhs ! (ikr,ikz)
  real(dp), DIMENSION(:,:), ALLOCATABLE :: laplace_sol ! (ikr,ikz)

  ! For FFT
  double complex, DIMENSION(:,:), ALLOCATABLE :: fft_filter ! (ikr,ikz)
  double complex, DIMENSION(:,:), ALLOCATABLE :: fft_input  ! (ikr,ikz)

  ! For WENO gradients
  real(dp), DIMENSION(:), ALLOCATABLE :: tmp1
  real(dp), DIMENSION(:), ALLOCATABLE :: tmp2
  real(dp), DIMENSION(:,:), ALLOCATABLE :: omegak_times_q3k


END MODULE array

