MODULE array

  !USE mumps_bsplines
  use prec_const
  implicit none

  ! Arrays to store the rhs, for time integration
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: moments_rhs ! (ipj,ikr,ikz,updatetlevel)

  ! Intermediate steps in rhs of equations
  !COMPLEX(dp), DIMENSION(:,:,:), ALLOCATABLE:: moments_Apl, moments_Bpl, moments_Cpl, moments_Dpl,&
  !                                        moments_Epl, moments_Fpl, moments_Gpl, moments_Hpl  
END MODULE array

