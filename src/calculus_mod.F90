MODULE calculus
  ! Routine to evaluate gradients, interpolation schemes and integrals
  USE basic
  USE prec_const
  USE grid
  USE parallel, ONLY: manual_0D_bcast
  IMPLICIT NONE
  REAL(dp), dimension(-2:2) :: dz_usu = &
   (/  onetwelfth, -twothird, &
                       0._dp, &
         twothird, -onetwelfth /) ! fd4 centered stencil
  REAL(dp), dimension(-2:1) :: dz_o2e = &
   (/ onetwentyfourth,-nineeighths, nineeighths,-onetwentyfourth /) ! fd4 odd to even stencil
  REAL(dp), dimension(-1:2) :: dz_e2o = &
   (/ onetwentyfourth,-nineeighths, nineeighths,-onetwentyfourth /) ! fd4 odd to even stencil
   REAL(dp), dimension(-2:2) :: dz2_usu = &
   (/-1.0_dp/12.0_dp, 4.0_dp/3.0_dp, -5.0_dp/2.0_dp, 4.0_dp/3.0_dp, -1.0_dp/12.0_dp /)! 2th derivative, 4th order (for parallel hypdiff)
   REAL(dp), dimension(-2:2) :: dz4_usu = &
   (/  1._dp, -4._dp, 6._dp, -4._dp, 1._dp /) ! 4th derivative, 2nd order (for parallel hypdiff)
  PUBLIC :: simpson_rule_z, interp_z, grad_z, grad_z4

CONTAINS

SUBROUTINE grad_z(target,f,ddzf)
  implicit none
  ! Compute the periodic boundary condition 4 points centered finite differences
  ! formula among staggered grid or not.
  ! not staggered : the derivative results must be on the same grid as the field
  !     staggered : the derivative is computed from a grid to the other
  INTEGER, INTENT(IN) :: target
  COMPLEX(dp),dimension(izgs:izge), intent(in)  :: f
  COMPLEX(dp),dimension( izs:ize ), intent(out) :: ddzf
  IF(Nz .GT. 3) THEN ! Cannot apply four points stencil on less than four points grid
    IF(SG) THEN
      IF(TARGET .EQ. 0) THEN
        CALL grad_z_o2e(f,ddzf)
      ELSE
        CALL grad_z_e2o(f,ddzf)
      ENDIF
    ELSE ! No staggered grid -> usual centered finite differences
      DO iz = izs,ize
       ddzf(iz) = dz_usu(-2)*f(iz-2) + dz_usu(-1)*f(iz-1) &
                 +dz_usu( 1)*f(iz+1) + dz_usu( 2)*f(iz+2)
      ENDDO
      ddzf(:) = ddzf(:) * inv_deltaz
    ENDIF
  ELSE
    ddzf = 0._dp
  ENDIF
CONTAINS
  SUBROUTINE grad_z_o2e(fo,ddzfe) ! Paruta 2018 eq (27)
    ! gives the gradient of a field from the odd grid to the even one
    implicit none
    COMPLEX(dp),dimension(izgs:izge), intent(in)  :: fo
    COMPLEX(dp),dimension( izs:ize ), intent(out) :: ddzfe !
    DO iz = izs,ize
     ddzfe(iz) = dz_o2e(-2)*fo(iz-2) + dz_o2e(-1)*fo(iz-1) &
                +dz_o2e( 0)*fo(iz  ) + dz_o2e( 1)*fo(iz+1)
    ENDDO
    ddzfe(:) = ddzfe(:) * inv_deltaz
  END SUBROUTINE grad_z_o2e

  SUBROUTINE grad_z_e2o(fe,ddzfo) ! n2v for Paruta 2018 eq (28)
    ! gives the gradient of a field from the even grid to the odd one
    implicit none
    COMPLEX(dp),dimension(izgs:izge), intent(in)  :: fe
    COMPLEX(dp),dimension( izs:ize ), intent(out) :: ddzfo
    DO iz = izs,ize
     ddzfo(iz) = dz_e2o(-1)*fe(iz-1) + dz_e2o(0)*fe(iz  ) &
                +dz_e2o( 1)*fe(iz+1) + dz_e2o(2)*fe(iz+2)
    ENDDO
    ddzfo(:) = ddzfo(:) * inv_deltaz
  END SUBROUTINE grad_z_e2o
END SUBROUTINE grad_z

SUBROUTINE grad_z2(f,ddz2f)
  implicit none
  ! Compute the second order fourth derivative for periodic boundary condition
  COMPLEX(dp),dimension(izgs:izge), intent(in)  :: f
  COMPLEX(dp),dimension( izs:ize ), intent(out) :: ddz2f
  IF(Nz .GT. 3) THEN ! Cannot apply four points stencil on less than four points grid
      DO iz = izs,ize
       ddz2f(iz) = dz2_usu(-2)*f(iz-2) + dz2_usu(-1)*f(iz-1) &
                  +dz2_usu( 0)*f(iz  )&
                  +dz2_usu( 1)*f(iz+1) + dz2_usu( 2)*f(iz+2)
      ENDDO
  ELSE
    ddz2f = 0._dp
  ENDIF
  ddz2f = ddz2f * inv_deltaz**2
END SUBROUTINE grad_z2


SUBROUTINE grad_z4(f,ddz4f)
  implicit none
  ! Compute the second order fourth derivative for periodic boundary condition
  COMPLEX(dp),dimension(izgs:izge), intent(in)  :: f
  COMPLEX(dp),dimension( izs:ize ), intent(out) :: ddz4f
  IF(Nz .GT. 3) THEN ! Cannot apply four points stencil on less than four points grid
      DO iz = izs,ize
       ddz4f(iz) = dz4_usu(-2)*f(iz-2) + dz4_usu(-1)*f(iz-1) &
                  +dz4_usu( 0)*f(iz)&
                  +dz4_usu( 1)*f(iz+1) + dz4_usu( 2)*f(iz+2)
      ENDDO
  ELSE
    ddz4f = 0._dp
  ENDIF
  ddz4f = ddz4f * inv_deltaz**4
END SUBROUTINE grad_z4


SUBROUTINE interp_z(target,f_in,f_out)
  ! Function meant to interpolate one field defined on a even/odd z into
  !  the other odd/even z grid.
  ! If Staggered Grid flag (SG) is false, returns identity
  implicit none
  INTEGER, intent(in) :: target ! target grid : 0 for even grid, 1 for odd
  COMPLEX(dp),dimension(izgs:izge), intent(in)  :: f_in
  COMPLEX(dp),dimension( izs:ize ), intent(out) :: f_out !
  IF(SG) THEN
    IF(target .EQ. 0) THEN
      CALL interp_o2e_z(f_in,f_out)
    ELSE
      CALL interp_e2o_z(f_in,f_out)
    ENDIF
  ELSE ! No staggered grid -> identity
    f_out(izs:ize) = f_in(izs:ize)
  ENDIF
CONTAINS
  SUBROUTINE interp_o2e_z(fo,fe)
   ! gives the value of a field from the odd grid to the even one
   implicit none
   COMPLEX(dp),dimension(izgs:izge), intent(in)  :: fo
   COMPLEX(dp),dimension( izs:ize ), intent(out) :: fe !
   ! 4th order interp
   DO iz = izs,ize
    fe(iz) = onesixteenth * (-fo(iz-2) + 9._dp*(fo(iz-1) + fo(iz)) - fo(iz+1))
   ENDDO
  END SUBROUTINE interp_o2e_z

  SUBROUTINE interp_e2o_z(fe,fo)
   ! gives the value of a field from the even grid to the odd one
   implicit none
   COMPLEX(dp),dimension(izgs:izge), intent(in)  :: fe
   COMPLEX(dp),dimension( izs:ize ), intent(out) :: fo
   ! 4th order interp
   DO iz = izs,ize
    fo(iz) = onesixteenth * (-fe(iz-1) + 9._dp*(fe(iz) + fe(iz+1)) - fe(iz+2))
   ENDDO
  END SUBROUTINE interp_e2o_z

END SUBROUTINE interp_z

SUBROUTINE simpson_rule_z(f,intf)
 ! integrate f(z) over z using the simpon's rule. Assume periodic boundary conditions (f(ize+1) = f(izs))
 !from molix BJ Frei
 implicit none
 complex(dp),dimension(izs:ize), intent(in) :: f
 COMPLEX(dp), intent(out) :: intf
 COMPLEX(dp)              :: buffer, local_int
 INTEGER :: root, i_

 IF(Nz .EQ. 1) THEN !2D zpinch simulations
   intf = f(izs)

 ELSE !3D fluxtube
   IF(mod(Nz,2) .ne. 0 ) THEN
      ERROR STOP '>> ERROR << Simpson rule: Nz must be an even number  !!!!'
   ENDIF
   ! Buil local sum using the weights of composite Simpson's rule
   local_int = 0._dp
   DO iz = izs,ize
      local_int = local_int + zweights_SR(iz)*f(iz)
   ENDDO
   buffer = local_int
   root = 0
   !Gather manually among the rank_z=0 processes and perform the sum
   intf = 0._dp
   IF (num_procs_z .GT. 1) THEN
       !! Everyone sends its local_sum to root = 0
       IF (rank_z .NE. root) THEN
           CALL MPI_SEND(buffer, 1 , MPI_DOUBLE_COMPLEX, root, 5678, comm_z, ierr)
       ELSE
           ! Recieve from all the other processes
           DO i_ = 0,num_procs_z-1
               IF (i_ .NE. rank_z) &
                   CALL MPI_RECV(buffer, 1 , MPI_DOUBLE_COMPLEX, i_, 5678, comm_z, MPI_STATUS_IGNORE, ierr)
                   intf = intf + buffer
           ENDDO
       ENDIF
       CALL manual_0D_bcast(intf)
   ELSE
     intf = local_int
   ENDIF
   intf = onethird*deltaz*intf
   ENDIF
END SUBROUTINE simpson_rule_z

END MODULE calculus
