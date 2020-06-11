module gradients
  use space_grid
  use prec_const
  implicit none

  ! Staggered grids are used, the v-grid of vpar is shifted by deltaz/2 compared to the normal n-grid of theta, temp, moments, phi.
  ! (z_i) on n-grid is (z_i+deltaz/2) of v-grid.
  ! When evaluating a gradient on another grid, we use formulas that interpolate directly.

  ! Gradients
  procedure(gradsub), pointer :: gradz => null()
  procedure(gradsub), pointer :: gradz_n2v => null()
  procedure(gradsub), pointer :: gradz_v2n => null()

  procedure(gradsub), pointer :: gradzm => null() ! For upwind scheme, other set of spatial derivatives, the "minus side"

  procedure(gradsub), pointer :: gradzz => null()


  ! Grid-grid interpolations
  procedure(gradsub), pointer :: interp_n2v => null()
  procedure(gradsub), pointer :: interp_v2n => null()



  !Interface for gradient subroutines 
  abstract interface
     subroutine gradsub( f , fp )
       use space_grid, only: izs,ize
       use prec_const
       implicit none

       real(dp), dimension(izs:ize), intent(in)  :: f
       real(dp), dimension(izs:ize), intent(out) :: fp

     end subroutine gradsub
  end interface


contains
  !To promote cleanliness, use .h include files
  !Please put things where they belong if you add new routines

  !Inside this file place formulas that are order independent,
  !e.g. sums of gradients and such

  !Note that eventually all the functions here should become private

#include "grad_fa2.h"
#include "grad_up2.h"
#include "grad_fd4.h"
#include "grad_we4.h"


  !In this subroutine, which should be called before any gradient is carried out,
  !we select the numerical scheme using function pointers
  subroutine set_gradient_scheme(gradient_scheme)

    use prec_const
    implicit none
    character(len=3),intent(in) :: gradient_scheme

    select case ( gradient_scheme )
! NOT IMPLEMENTED YET, need finite difference order 2 formulas for staggered grids
!    case ( 'fa2' )
!       gradz => gradz_fa2
!       gradz_n2v => gradz_n2v_fa2
!       gradz_v2n => gradz_v2n_fa2
!       gradzz => gradzz_fa2
!    case ( 'up2' )
!       gradz  => gradz_up2_plus
!       gradzm  => gradz_up2_minus
    case ('fd4')
       ! Gradients
!       gradz => gradz_forfd4
       gradz => gradz_fd4
       gradz_n2v => gradz_n2v_fd4
       gradz_v2n => gradz_v2n_fd4
       gradzz => gradzz_fd4
       
       ! Interpolation
       interp_n2v => interp_n2v_fd4
       interp_v2n => interp_v2n_fd4

    case ('we4') ! WENO weighted essentially non-oscillatory order 4, meaning r=3
       ! Gradients
       gradz => gradz_we4
       gradz_n2v => gradz_n2v_fd4
       gradz_v2n => gradz_v2n_fd4
       gradzz => gradzz_fd4

       ! Interpolation
       interp_n2v => interp_n2v_fd4
       interp_v2n => interp_v2n_fd4
	

    case default
       write (*,*) "Unknown numerical scheme ", gradient_scheme
    end select

  end subroutine set_gradient_scheme

end module gradients
