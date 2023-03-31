SUBROUTINE compute_gradients
    !   This routine compute the gradients of the moments without copying or slices.
    ! It should be faster than using a routine taking a slice as argument (see 
    ! calculus_mod) since it avoid copying.
    USE fields, ONLY: moments
    USE grid,   ONLY: local_nz, ngz, inv_deltaz
    USE prec_const, ONLY: dp
    IMPLICIT NONE

    REAL(dp), dimension(-2:2) :: dz_usu = &
     (/  1._dp/12._dp, -2._dp/3._dp, 0._dp, 2._dp/3._dp, -1._dp/12._dp /) ! fd4 centered stencil
    REAL(dp), dimension(-2:1) :: dz_o2e = &
     (/ 1._dp/24._dp,-9._dp/8._dp, 9._dp/8._dp,-1._dp/24._dp /) ! fd4 odd to even stencil
    REAL(dp), dimension(-1:2) :: dz_e2o = &
     (/ 1._dp/24._dp,-9._dp/8._dp, 9._dp/8._dp,-1._dp/24._dp /) ! fd4 odd to even stencil
     REAL(dp), dimension(-2:2) :: dz2_usu = &
     (/-1._dp/12._dp, 4._dp/3._dp, -5._dp/2._dp, 4._dp/3._dp, -1._dp/12._dp /)! 2th derivative, 4th order (for parallel hypdiff)
     REAL(dp), dimension(-2:2) :: dz4_usu = &
     (/  1._dp, -4._dp, 6._dp, -4._dp, 1._dp /) ! 4th derivative, 2nd order (for parallel hypdiff)
     REAL(dp), dimension(-2:1) :: iz_o2e = &
     (/ -1._dp/16._dp, 9._dp/16._dp, 9._dp/16._dp, -1._dp/16._dp /) ! grid interpolation, 4th order, odd to even
     REAL(dp), dimension(-1:2) :: iz_e2o = &
     (/ -1._dp/16._dp, 9._dp/16._dp, 9._dp/16._dp, -1._dp/16._dp /) ! grid interpolation, 4th order, even to odd
    PUBLIC :: simpson_rule_z, interp_z, grad_z, grad_z4
  
IMPLICIT NONE

END SUBROUTINE compute_gradients