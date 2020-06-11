SUBROUTINE gradz_fd4( f , f_z ) ! implemented for periodic boundary conditions

  use prec_const
  use space_grid
  implicit none
  real(dp), dimension(izs:ize),intent(in) :: f
  real(dp), dimension(izs:ize),intent(out) :: f_z
  real(dp), dimension(-2:2) :: coef ! fd4 stencil
  integer :: iz

  coef = (/ deltazi*onetwelfth, &
            -deltazi*2._dp*onethird, &
            0._dp, &
            deltazi*2._dp*onethird, &
            -deltazi*onetwelfth /) 

  f_z(izs)   = coef(-2)*f(ize-1) + coef(-1)*f(ize) + coef(1)*f(izs+1) + coef(2)*f(izs+2)
  f_z(izs+1) = coef(-2)*f(ize) + coef(-1)*f(izs) + coef(1)*f(izs+2) + coef(2)*f(izs+3)
  do iz=izs+2,ize-2
    f_z(iz)  = coef(-2)*f(iz-2) + coef(-1)*f(iz-1) + coef(1)*f(iz+1) + coef(2)*f(iz+2)  
  end do
  f_z(ize-1) = coef(-2)*f(ize-3) + coef(-1)*f(ize-2) + coef(1)*f(ize) + coef(2)*f(izs)
  f_z(ize)   = coef(-2)*f(ize-2) + coef(-1)*f(ize-1) + coef(1)*f(izs) + coef(2)*f(izs+1)

END SUBROUTINE gradz_fd4


SUBROUTINE gradz_n2v_fd4( f , f_z ) ! implemented for periodic boundary conditions

  use prec_const
  use space_grid
  implicit none
  real(dp), dimension(izs:ize),intent(in) :: f
  real(dp), dimension(izs:ize),intent(out) :: f_z
  real(dp), dimension(-2:1) :: coef ! fd4 stencil
  integer :: iz

  coef = (/ deltazi*onetwentyfourth, &
            -deltazi*nineeighths, &
            deltazi*nineeighths, &
            -deltazi*onetwentyfourth /) 

  f_z(izs)   = coef(-2)*f(ize-1) + coef(-1)*f(ize) + coef(0)*f(izs) + coef(1)*f(izs+1)
  f_z(izs+1) = coef(-2)*f(ize) + coef(-1)*f(izs) + coef(0)*f(izs+1) + coef(1)*f(izs+2)
  do iz=izs+2,ize-1
    f_z(iz)  = coef(-2)*f(iz-2) + coef(-1)*f(iz-1) + coef(0)*f(iz) + coef(1)*f(iz+1)  
  end do
  f_z(ize)   = coef(-2)*f(ize-2) + coef(-1)*f(ize-1) + coef(0)*f(ize) + coef(1)*f(izs)

END SUBROUTINE gradz_n2v_fd4


SUBROUTINE gradz_v2n_fd4( f , f_z ) ! implemented for periodic boundary conditions

  use prec_const
  use space_grid
  implicit none
  real(dp), dimension(izs:ize),intent(in) :: f
  real(dp), dimension(izs:ize),intent(out) :: f_z
  real(dp), dimension(-1:2) :: coef ! fd4 stencil
  integer :: iz

  coef = (/ deltazi*onetwentyfourth, &
            -deltazi*nineeighths, &
            deltazi*nineeighths, &
            -deltazi*onetwentyfourth /) 

  f_z(izs) = coef(-1)*f(ize) + coef(0)*f(izs) + coef(1)*f(izs+1) + coef(2)*f(izs+2)
  do iz=izs+1,ize-2
    f_z(iz) = coef(-1)*f(iz-1) + coef(0)*f(iz) + coef(1)*f(iz+1) + coef(2)*f(iz+2)
  end do
  f_z(ize-1) = coef(-1)*f(ize-2) + coef(0)*f(ize-1) + coef(1)*f(ize) + coef(2)*f(izs)
  f_z(ize)   = coef(-1)*f(ize-1) + coef(0)*f(ize) + coef(1)*f(izs) + coef(2)*f(izs+1)

END SUBROUTINE gradz_v2n_fd4



SUBROUTINE gradzz_fd4( f , f_z ) ! implemented for periodic boundary conditions

  use prec_const
  use space_grid
  implicit none
  real(dp), dimension(izs:ize),intent(in) :: f
  real(dp), dimension(izs:ize),intent(out) :: f_z
  real(dp), dimension(-2:2) :: coef ! fd4 stencil
  integer :: iz

  coef = (/ -deltazisq*onetwelfth, &
            deltazisq*4._dp*onethird, &
            -deltazisq*2.5_dp, &
            deltazisq*4._dp*onethird, &
            -deltazisq*onetwelfth /) 

  f_z(izs) = coef(-2)*f(ize-1) + coef(-1)*f(ize) + coef(0)*f(izs) + coef(1)*f(izs+1) + coef(2)*f(izs+2)
  f_z(izs+1) = coef(-2)*f(ize) + coef(-1)*f(izs) + coef(0)*f(izs+1) + coef(1)*f(izs+2) + coef(2)*f(izs+3)
  do iz=izs+2,ize-2
    f_z(iz) = coef(-2)*f(iz-2) + coef(-1)*f(iz-1) + coef(0)*f(iz) + coef(1)*f(iz+1) + coef(2)*f(iz+2)  
  end do
  f_z(ize-1) = coef(-2)*f(ize-3) + coef(-1)*f(ize-2) + coef(0)*f(ize-1) + coef(1)*f(ize) + coef(2)*f(izs)
  f_z(ize)   = coef(-2)*f(ize-2) + coef(-1)*f(ize-1) + coef(0)*f(ize) + coef(1)*f(izs) + coef(2)*f(izs+1)

END SUBROUTINE gradzz_fd4



! Interpolation

SUBROUTINE interp_n2v_fd4( f_n , f_v ) ! implemented for periodic boundary conditions

  use prec_const
  use space_grid
  implicit none
  real(dp), dimension(izs:ize),intent(in) :: f_n
  real(dp), dimension(izs:ize),intent(out) :: f_v
  real(dp), dimension(-2:1) :: coef ! fd4 stencil
  integer :: iz

  coef = (/ -onesixteenth, &
            ninesixteenths, &
            ninesixteenths, &
            -onesixteenth /) 

  f_v(izs)   = coef(-2)*f_n(ize-1) + coef(-1)*f_n(ize) + coef(0)*f_n(izs) + coef(1)*f_n(izs+1)
  f_v(izs+1) = coef(-2)*f_n(ize) + coef(-1)*f_n(izs) + coef(0)*f_n(izs+1) + coef(1)*f_n(izs+2)
  do iz=izs+2,ize-1
    f_v(iz)  = coef(-2)*f_n(iz-2) + coef(-1)*f_n(iz-1) + coef(0)*f_n(iz) + coef(1)*f_n(iz+1)  
  end do
  f_v(ize)   = coef(-2)*f_n(ize-2) + coef(-1)*f_n(ize-1) + coef(0)*f_n(ize) + coef(1)*f_n(izs)

END SUBROUTINE interp_n2v_fd4



SUBROUTINE interp_v2n_fd4( f_v , f_n ) ! implemented for periodic boundary conditions

  use prec_const
  use space_grid
  implicit none
  real(dp), dimension(izs:ize),intent(in) :: f_v
  real(dp), dimension(izs:ize),intent(out) :: f_n
  real(dp), dimension(-1:2) :: coef ! fd4 stencil
  integer :: iz

  coef = (/ -onesixteenth, &
            ninesixteenths, &
            ninesixteenths, &
            -onesixteenth /) 

  f_n(izs) = coef(-1)*f_v(ize) + coef(0)*f_v(izs) + coef(1)*f_v(izs+1) + coef(2)*f_v(izs+2)
  do iz=izs+1,ize-2
    f_n(iz) = coef(-1)*f_v(iz-1) + coef(0)*f_v(iz) + coef(1)*f_v(iz+1) + coef(2)*f_v(iz+2)
  end do
  f_n(ize-1) = coef(-1)*f_v(ize-2) + coef(0)*f_v(ize-1) + coef(1)*f_v(ize) + coef(2)*f_v(izs)
  f_n(ize)   = coef(-1)*f_v(ize-1) + coef(0)*f_v(ize) + coef(1)*f_v(izs) + coef(2)*f_v(izs+1)

END SUBROUTINE interp_v2n_fd4


