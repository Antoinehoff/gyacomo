SUBROUTINE gradz_fa2( f , f_z )

  use prec_const
  implicit none
  real(dp), dimension(izs:ize),intent(in) :: f
  real(dp), dimension(izs:ize),intent(out) :: f_z
  integer :: iz

  f_z(izs) = deltazih * ( f(izs+1) - f(ize) )
  do iz=izs+1,ize-1
    f_z(iz) = deltazih * ( f(iz+1) - f(iz-1) )
  end do
  f_z(ize) = deltazih * ( f(izs) - f(ize-1) )

END SUBROUTINE gradz_fa2


SUBROUTINE gradzz_fa2( f , f_z )

  use prec_const
  implicit none
  real(dp), dimension(izs:ize),intent(in) :: f
  real(dp), dimension(izs:ize),intent(out) :: f_z
  integer :: iz

  f_z(izs) = deltazisq * ( f(izs+1) - 2._dp*f(izs) + f(ize) )
  do iz=izs+1,ize-1
    f_z(iz) = deltazisq * ( f(iz+1) - 2._dp*f(iz) + f(iz-1) )
  end do
  f_z(ize) = deltazisq * ( f(izs) - 2._dp*f(ize) + f(ize-1) )

END SUBROUTINE gradzz_fa2
