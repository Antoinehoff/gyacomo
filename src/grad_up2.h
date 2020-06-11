! For the upwind scheme, the gradients using the two possible sides 
SUBROUTINE gradz_up2_plus( f , f_z )

  use prec_const
  implicit none
  real(dp), dimension(izs:ize),intent(in) :: f
  real(dp), dimension(izs:ize),intent(out) :: f_z
  integer :: iz

  do iz=izs,ize-2
    f_z(iz) = deltazih * ( -3._dp*f(iz) + 4._dp*f(iz+1) -f(iz+2) )
  end do
  f_z(ize-1) = deltazih * ( -3._dp*f(ize-1) + 4._dp*f(ize) -f(izs) )
  f_z(ize) = deltazih * ( -3._dp*f(ize) + 4._dp*f(izs) -f(izs+1) )

END SUBROUTINE gradz_up2_plus


SUBROUTINE gradz_up2_minus( f , f_z )

  use prec_const
  implicit none
  real(dp), dimension(izs:ize),intent(in) :: f
  real(dp), dimension(izs:ize),intent(out) :: f_z
  integer :: iz

  f_z(izs) = deltazih * ( f(ize-1) - 4._dp*f(ize) + 3._dp*f(izs) )
  f_z(izs+1) = deltazih * ( f(ize) - 4._dp*f(izs) + 3._dp*f(izs+1) )
  do iz=izs+2,ize
    f_z(iz) = deltazih * ( f(iz-2) - 4._dp*f(iz-1) + 3._dp*f(iz) )
  end do

END SUBROUTINE gradz_up2_minus
