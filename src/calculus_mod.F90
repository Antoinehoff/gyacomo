MODULE calculus
  ! Routine to evaluate gradients, interpolation schemes and integrals
  USE basic
  USE prec_const
  USE grid
  IMPLICIT NONE
  REAL(dp), dimension(-2:2) :: dz_usu = &
   (/  onetwelfth, -twothird, &
                       0._dp, &
         twothird, -onetwelfth /) ! fd4 centered stencil
  REAL(dp), dimension(-2:1) :: dz_o2e = &
   (/ onetwentyfourth,-nineeighths, nineeighths,-onetwentyfourth /) ! fd4 odd to even stencil
  REAL(dp), dimension(-1:2) :: dz_e2o = &
   (/ onetwentyfourth,-nineeighths, nineeighths,-onetwentyfourth /) ! fd4 odd to even stencil

  PUBLIC :: simpson_rule_z, interp_z, grad_z

CONTAINS

SUBROUTINE grad_z(target,f,ddzf)
  implicit none
  ! Compute the periodic boundary condition 4 points centered finite differences
  ! formula among staggered grid or not.
  ! not staggered : the derivative results must be on the same grid as the field
  !     staggered : the derivative is computed from a grid to the other
  INTEGER, INTENT(IN) :: target
  COMPLEX(dp),dimension( izs:ize ), intent(in)  :: f
  COMPLEX(dp),dimension( izs:ize ), intent(out) :: ddzf
  IF(Nz .GT. 3) THEN ! Cannot apply four points stencil on less than four points grid
    IF(SG) THEN
      IF(TARGET .EQ. 0) THEN
        CALL grad_z_o2e(f,ddzf)
      ELSE
        CALL grad_z_e2o(f,ddzf)
      ENDIF
    ELSE ! No staggered grid -> usual centered finite differences
      DO iz = 3,Nz-2
       ddzf(iz) = dz_usu(-2)*f(iz-2) + dz_usu(-1)*f(iz-1) &
                 +dz_usu( 1)*f(iz+1) + dz_usu( 2)*f(iz+2)
      ENDDO
      ! Periodic boundary conditions
      ddzf(Nz-1) = dz_usu(-2)*f(Nz-3) + dz_usu(-1)*f(Nz-2) &
                  +dz_usu(+1)*f(Nz  ) + dz_usu(+2)*f(1   )
      ddzf(Nz  ) = dz_usu(-2)*f(Nz-2) + dz_usu(-1)*f(Nz-1) &
                  +dz_usu(+1)*f(1   ) + dz_usu(+2)*f(2   )
      ddzf(1   ) = dz_usu(-2)*f(Nz-1) + dz_usu(-1)*f(Nz  ) &
                  +dz_usu(+1)*f(2   ) + dz_usu(+2)*f(3   )
      ddzf(2   ) = dz_usu(-2)*f(Nz  ) + dz_usu(-1)*f(1   ) &
                  +dz_usu(+1)*f(3   ) + dz_usu(+2)*f(4)
    ENDIF
  ELSE
    ddzf = 0._dp
  ENDIF
END SUBROUTINE grad_z

SUBROUTINE grad_z_o2e(fo,dzfe) ! Paruta 2018 eq (27)
  ! gives the value of a field from the odd grid to the even one
  implicit none
  COMPLEX(dp),dimension(1:Nz), intent(in)  :: fo
  COMPLEX(dp),dimension(1:Nz), intent(out) :: dzfe !
  DO iz = 3,Nz-1
   dzfe(iz) = dz_o2e(-2)*fo(iz-2) + dz_o2e(-1)*fo(iz-1) &
             +dz_o2e( 0)*fo(iz  ) + dz_o2e( 1)*fo(iz+1)
  ENDDO
  ! Periodic boundary conditions
  dzfe(Nz) = dz_o2e(-2)*fo(Nz-2) + dz_o2e(-1)*fo(Nz-1) &
            +dz_o2e( 0)*fo(Nz  ) + dz_o2e( 1)*fo(1)
  dzfe(1 ) = dz_o2e(-2)*fo(Nz-1) + dz_o2e(-1)*fo(Nz  ) &
            +dz_o2e( 0)*fo(1   ) + dz_o2e( 1)*fo(2)
  dzfe(2 ) = dz_o2e(-2)*fo(Nz  ) + dz_o2e(-1)*fo(1   ) &
            +dz_o2e( 0)*fo(2   ) + dz_o2e( 1)*fo(3)
END SUBROUTINE grad_z_o2e

SUBROUTINE grad_z_e2o(fe,dzfo)
  ! gives the value of a field from the even grid to the odd one
  implicit none
  COMPLEX(dp),dimension(1:Nz), intent(in)  :: fe
  COMPLEX(dp),dimension(1:Nz), intent(out) :: dzfo
  DO iz = 2,Nz-2
   dzfo(iz) = dz_e2o(-1)*fe(iz-1) + dz_e2o(0)*fe(iz  ) &
             +dz_e2o( 1)*fe(iz+1) + dz_e2o(2)*fe(iz+2)
  ENDDO
  ! Periodic boundary conditions
  dzfo(Nz-1) = dz_e2o(-1)*fe(Nz-2) + dz_e2o(0)*fe(Nz-1) &
              +dz_e2o( 1)*fe(Nz  ) + dz_e2o(2)*fe(1   )
  dzfo(Nz  ) = dz_e2o(-1)*fe(Nz-1) + dz_e2o(0)*fe(Nz  ) &
              +dz_e2o( 1)*fe(1   ) + dz_e2o(2)*fe(2   )
  dzfo(1   ) = dz_e2o(-1)*fe(Nz  ) + dz_e2o(0)*fe(1   ) &
              +dz_e2o( 1)*fe(2   ) + dz_e2o(2)*fe(3   )
END SUBROUTINE grad_z_e2o


SUBROUTINE interp_z(target,f_in,f_out)
  ! Function meant to interpolate one field defined on a even/odd z into
  !  the other odd/even z grid.
  ! If Staggered Grid flag (SG) is false, returns identity
  implicit none
  INTEGER, intent(in) :: target ! target grid : 0 for even grid, 1 for odd
  COMPLEX(dp),dimension(1:Nz), intent(in)  :: f_in
  COMPLEX(dp),dimension(1:Nz), intent(out) :: f_out !
  IF(SG) THEN
    IF(TARGET .EQ. 0) THEN
      CALL interp_o2e_z(f_in,f_out)
    ELSE
      CALL interp_e2o_z(f_in,f_out)
    ENDIF
  ELSE ! No staggered grid -> identity
    f_out(:) = f_in(:)
  ENDIF
END SUBROUTINE interp_z

SUBROUTINE interp_o2e_z(fo,fe)
 ! gives the value of a field from the odd grid to the even one
 implicit none
 COMPLEX(dp),dimension(1:Nz), intent(in)  :: fo
 COMPLEX(dp),dimension(1:Nz), intent(out) :: fe !
 DO iz = 2,Nz
   fe(iz) = 0.5_dp*(fo(iz)+fo(iz-1))
 ENDDO
 ! Periodic boundary conditions
 fe(1) = 0.5_dp*(fo(1) + fo(Nz))
END SUBROUTINE interp_o2e_z

SUBROUTINE interp_e2o_z(fe,fo)
 ! gives the value of a field from the even grid to the odd one
 implicit none
 COMPLEX(dp),dimension(1:Nz), intent(in)  :: fe
 COMPLEX(dp),dimension(1:Nz), intent(out) :: fo
 DO iz = 1,Nz-1
   fo(iz) = 0.5_dp*(fe(iz+1)+fe(iz))
 ENDDO
 ! Periodic boundary conditions
 fo(Nz) = 0.5_dp*(fe(1) + fe(Nz))
END SUBROUTINE interp_e2o_z


SUBROUTINE simpson_rule_z(f,intf)
 ! integrate f(z) over z using the simpon's rule. Assume periodic boundary conditions (f(ize+1) = f(izs))
 !from molix BJ Frei
 implicit none
 complex(dp),dimension(izs:ize), intent(in) :: f
 COMPLEX(dp), intent(out) :: intf
 COMPLEX(dp) :: buff_
 IF(Nz .GT. 1) THEN
   IF(mod(Nz,2) .ne. 0 ) THEN
      ERROR STOP 'Simpson rule: Nz must be an even number  !!!!'
   ENDIF
   buff_ = 0._dp
   DO iz = izs, Nz/2
      IF(iz .eq. Nz/2) THEN ! ... iz = ize
         buff_ = buff_ + (f(izs) + 4._dp*f(ize) + f(ize-1 ))
      ELSE
         buff_ = buff_ + (f(2*iz+1) + 4._dp*f(2*iz) + f(2*iz-1 ))
      ENDIF
   ENDDO
   intf = buff_*deltaz/3._dp
 ELSE
   intf = f(izs)
 ENDIF
END SUBROUTINE simpson_rule_z

END MODULE calculus
