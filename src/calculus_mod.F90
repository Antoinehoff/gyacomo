MODULE calculus
  ! Routine to evaluate gradients, interpolation schemes and integrals
  USE prec_const, ONLY: xp
  IMPLICIT NONE
  REAL(xp), dimension(-2:2) :: dz_usu = &
   (/  1._xp/12._xp, -2._xp/3._xp, 0._xp, 2._xp/3._xp, -1._xp/12._xp /) ! fd4 centered stencil
  REAL(xp), dimension(-2:1) :: dz_o2e = &
   (/ 1._xp/24._xp,-9._xp/8._xp, 9._xp/8._xp,-1._xp/24._xp /) ! fd4 odd to even stencil
  REAL(xp), dimension(-1:2) :: dz_e2o = &
   (/ 1._xp/24._xp,-9._xp/8._xp, 9._xp/8._xp,-1._xp/24._xp /) ! fd4 odd to even stencil
   REAL(xp), dimension(-2:2) :: dz2_usu = &
   (/-1._xp/12._xp, 4._xp/3._xp, -5._xp/2._xp, 4._xp/3._xp, -1._xp/12._xp /)! 2th derivative, 4th order (for parallel hypdiff)
   REAL(xp), dimension(-2:2) :: dz4_usu = &
   (/  1._xp, -4._xp, 6._xp, -4._xp, 1._xp /) ! 4th derivative, 2nd order (for parallel hypdiff)
   REAL(xp), dimension(-2:1) :: iz_o2e = &
   (/ -1._xp/16._xp, 9._xp/16._xp, 9._xp/16._xp, -1._xp/16._xp /) ! grid interpolation, 4th order, odd to even
   REAL(xp), dimension(-1:2) :: iz_e2o = &
   (/ -1._xp/16._xp, 9._xp/16._xp, 9._xp/16._xp, -1._xp/16._xp /) ! grid interpolation, 4th order, even to odd
  PUBLIC :: simpson_rule_z, interp_z, grad_z, grad_z4, grad_z_5D, grad_z4_5D

CONTAINS

SUBROUTINE grad_z(target,local_nz,ngz,inv_deltaz,f,ddzf)
  implicit none
  ! Compute the periodic boundary condition 4 points centered finite differences
  ! formula among staggered grid or not.
  ! not staggered : the derivative results must be on the same grid as the field
  !     staggered : the derivative is computed from a grid to the other
  INTEGER,  INTENT(IN) :: target, local_nz, ngz
  REAL(xp), INTENT(IN) :: inv_deltaz
  COMPLEX(xp),dimension(local_nz+ngz), INTENT(IN)  :: f
  COMPLEX(xp),dimension(local_nz),     INTENT(OUT) :: ddzf
  INTEGER :: iz
  IF(ngz .GT. 3) THEN ! Cannot apply four points stencil on less than four points grid
    SELECT CASE(TARGET)
    CASE(1)
      CALL grad_z_o2e(local_nz,ngz,inv_deltaz,f,ddzf)
    CASE(2)
      CALL grad_z_e2o(local_nz,ngz,inv_deltaz,f,ddzf)
    CASE DEFAULT ! No staggered grid -> usual centered finite differences
      DO iz = 1,local_nz
       ddzf(iz) = dz_usu(-2)*f(iz  ) + dz_usu(-1)*f(iz+1) &
                 +dz_usu( 0)*f(iz+2) &
                 +dz_usu( 1)*f(iz+3) + dz_usu( 2)*f(iz+4)
      ENDDO
      ddzf(:) = ddzf(:) * inv_deltaz
    END SELECT
  ELSE
    ddzf = 0._xp
  ENDIF
CONTAINS
  SUBROUTINE grad_z_o2e(local_nz,ngz,inv_deltaz,fo,ddzfe) ! Paruta 2018 eq (27)
    ! gives the gradient of a field from the odd grid to the even one
    implicit none
    INTEGER,  INTENT(IN) :: local_nz, ngz
    REAL(xp), INTENT(IN) :: inv_deltaz
    COMPLEX(xp),dimension(local_nz+ngz), INTENT(IN)  :: fo
    COMPLEX(xp),dimension(local_nz),     INTENT(OUT) :: ddzfe !
    INTEGER :: iz
    DO iz = 1,local_nz
     ddzfe(iz) = dz_o2e(-2)*fo(iz  ) + dz_o2e(-1)*fo(iz+1) &
                +dz_o2e( 0)*fo(iz+2) + dz_o2e( 1)*fo(iz+3)
    ENDDO
    ddzfe(:) = ddzfe(:) * inv_deltaz
  END SUBROUTINE grad_z_o2e

  SUBROUTINE grad_z_e2o(local_nz,ngz,inv_deltaz,fe,ddzfo) ! n2v for Paruta 2018 eq (28)
    ! gives the gradient of a field from the even grid to the odd one
    implicit none
    INTEGER,  INTENT(IN) :: local_nz, ngz
    REAL(xp), INTENT(IN) :: inv_deltaz
    COMPLEX(xp),dimension(local_nz+ngz), INTENT(IN)  :: fe
    COMPLEX(xp),dimension(local_nz),     INTENT(OUT) :: ddzfo
    INTEGER :: iz
    DO iz = 1,local_nz
     ddzfo(iz) = dz_e2o(-1)*fe(iz+1) + dz_e2o(0)*fe(iz+2) &
                +dz_e2o( 1)*fe(iz+3) + dz_e2o(2)*fe(iz+4)
    ENDDO
    ddzfo(:) = ddzfo(:) * inv_deltaz
  END SUBROUTINE grad_z_e2o
END SUBROUTINE grad_z

SUBROUTINE grad_z_5D(local_nz,ngz,inv_deltaz,f,ddzf)
  implicit none
  ! Compute the periodic boundary condition 4 points centered finite differences
  ! formula among staggered grid or not.
  ! not staggered : the derivative results must be on the same grid as the field
  !     staggered : the derivative is computed from a grid to the other
  INTEGER,  INTENT(IN) :: local_nz, ngz
  REAL(xp), INTENT(IN) :: inv_deltaz
  COMPLEX(xp),dimension(:,:,:,:,:,:), INTENT(IN)  :: f
  COMPLEX(xp),dimension(:,:,:,:,:,:),     INTENT(OUT) :: ddzf
  INTEGER :: iz
  IF(ngz .GT. 3) THEN ! Cannot apply four points stencil on less than four points grid
    DO iz = 1,local_nz
      ddzf(:,:,:,:,:,iz) = inv_deltaz*&
        (dz_usu(-2)*f(:,:,:,:,:,iz  ) + dz_usu(-1)*f(:,:,:,:,:,iz+1) &
        +dz_usu( 0)*f(:,:,:,:,:,iz+2) &
        +dz_usu( 1)*f(:,:,:,:,:,iz+3) + dz_usu( 2)*f(:,:,:,:,:,iz+4))
    ENDDO
  ELSE
    ddzf = 0._xp
  ENDIF
END SUBROUTINE grad_z_5D


SUBROUTINE grad_z2(local_nz,ngz,inv_deltaz,f,ddz2f)
  ! Compute the second order fourth derivative for periodic boundary condition
  implicit none
  INTEGER, INTENT(IN)  :: local_nz, ngz
  REAL(xp), INTENT(IN) :: inv_deltaz
  COMPLEX(xp),dimension(local_nz+ngz), INTENT(IN)  :: f
  COMPLEX(xp),dimension(local_nz),     INTENT(OUT) :: ddz2f
  INTEGER :: iz
  IF(ngz .GT. 3) THEN ! Cannot apply four points stencil on less than four points grid
      DO iz = 1,local_nz
       ddz2f(iz) = dz2_usu(-2)*f(iz  ) + dz2_usu(-1)*f(iz+1) &
                  +dz2_usu( 0)*f(iz+2)&
                  +dz2_usu( 1)*f(iz+3) + dz2_usu( 2)*f(iz+4)
      ENDDO
  ELSE
    ddz2f = 0._xp
  ENDIF
  ddz2f = ddz2f * inv_deltaz**2
END SUBROUTINE grad_z2

SUBROUTINE grad_z4_5D(local_nz,ngz,inv_deltaz,f,ddz4f)
  ! Compute the second order fourth derivative for periodic boundary condition
  implicit none
  INTEGER,  INTENT(IN) :: local_nz, ngz
  REAL(xp), INTENT(IN) :: inv_deltaz
  COMPLEX(xp),dimension(:,:,:,:,:,:), INTENT(IN)  :: f
  COMPLEX(xp),dimension(:,:,:,:,:,:),     INTENT(OUT) :: ddz4f
  INTEGER :: iz
  IF(ngz .GT. 3) THEN ! Cannot apply four points stencil on less than four points grid
      DO iz = 1,local_nz
       ddz4f(:,:,:,:,:,iz) = inv_deltaz**4*&
          (dz4_usu(-2)*f(:,:,:,:,:,iz  ) + dz4_usu(-1)*f(:,:,:,:,:,iz+1) &
          +dz4_usu( 0)*f(:,:,:,:,:,iz+2)&
          +dz4_usu( 1)*f(:,:,:,:,:,iz+3) + dz4_usu( 2)*f(:,:,:,:,:,iz+4))
      ENDDO
  ELSE
    ddz4f = 0._xp
  ENDIF
END SUBROUTINE grad_z4_5D

SUBROUTINE grad_z4(local_nz,ngz,inv_deltaz,f,ddz4f)
  ! Compute the second order fourth derivative for periodic boundary condition
  implicit none
  INTEGER,  INTENT(IN) :: local_nz, ngz
  REAL(xp), INTENT(IN) :: inv_deltaz
  COMPLEX(xp),dimension(local_nz+ngz), INTENT(IN)  :: f
  COMPLEX(xp),dimension(local_nz),     INTENT(OUT) :: ddz4f
  INTEGER :: iz
  IF(ngz .GT. 3) THEN ! Cannot apply four points stencil on less than four points grid
      DO iz = 1,local_nz
       ddz4f(iz) = dz4_usu(-2)*f(iz  ) + dz4_usu(-1)*f(iz+1) &
                  +dz4_usu( 0)*f(iz+2)&
                  +dz4_usu( 1)*f(iz+3) + dz4_usu( 2)*f(iz+4)
      ENDDO
  ELSE
    ddz4f = 0._xp
  ENDIF
  ddz4f = ddz4f * inv_deltaz**4
END SUBROUTINE grad_z4


SUBROUTINE interp_z(target,local_nz,ngz,f_in,f_out)
  ! Function meant to interpolate one field defined on a even/odd z into
  !  the other odd/even z grid.
  ! If Staggered Grid flag (SG) is false, returns identity
  implicit none
  INTEGER, INTENT(IN) :: local_nz, ngz
  INTEGER, intent(in) :: target ! target grid : 0 for even grid, 1 for odd
  COMPLEX(xp),dimension(local_nz+ngz), INTENT(IN)  :: f_in
  COMPLEX(xp),dimension(local_nz),     INTENT(OUT) :: f_out
  SELECT CASE(TARGET)
  CASE(1) ! output on even grid
    CALL interp_o2e_z(local_nz,ngz,f_in,f_out)
  CASE(2) ! output on odd grid
    CALL interp_e2o_z(local_nz,ngz,f_in,f_out)
  CASE DEFAULT ! No staggered grid -> usual centered finite differences
    f_out = f_in((1+ngz/2):(local_nz+ngz/2))
  END SELECT
CONTAINS
  SUBROUTINE interp_o2e_z(local_nz, ngz,fo,fe)
   ! gives the value of a field from the odd grid to the even one
   implicit none
   INTEGER, INTENT(IN) :: local_nz, ngz
   COMPLEX(xp),dimension(local_nz+ngz), INTENT(IN)  :: fo
   COMPLEX(xp),dimension(local_nz),     INTENT(OUT) :: fe
   INTEGER :: iz
   ! 4th order interp
   DO iz = 1,local_nz
     fe(iz) = iz_o2e(-2)*fo(iz )  + iz_o2e(-1)*fo(iz+1) &
            + iz_o2e( 0)*fo(iz+2) + iz_o2e( 1)*fo(iz+3)
   ENDDO
  END SUBROUTINE interp_o2e_z

  SUBROUTINE interp_e2o_z(local_nz, ngz,fe,fo)
   ! gives the value of a field from the even grid to the odd one
   implicit none
   INTEGER, INTENT(IN) :: local_nz, ngz
   COMPLEX(xp),dimension(local_nz+ngz), INTENT(IN)  :: fe
   COMPLEX(xp),dimension(local_nz),     INTENT(OUT) :: fo
   INTEGER :: iz
   ! 4th order interp
   DO iz = 1,local_nz
     fo(iz) = iz_e2o(-1)*fe(iz+1) + iz_e2o( 0)*fe(iz+2) &
            + iz_e2o( 1)*fe(iz+3) + iz_e2o( 2)*fe(iz+4)
   ENDDO
  END SUBROUTINE interp_e2o_z

END SUBROUTINE interp_z

SUBROUTINE simpson_rule_z(local_nz,zweights_SR,f,intf)
  ! integrate f(z) over z using the simpon's rule. Assume periodic boundary conditions (f(ize+1) = f(izs))
  !from molix BJ Frei
  USE prec_const, ONLY: xp, onethird, mpi_xp_c
  USE parallel,   ONLY: num_procs_z, rank_z, comm_z, manual_0D_bcast
  USE mpi
  implicit none
  INTEGER, INTENT(IN) :: local_nz
  REAL(xp),   dimension(local_nz), intent(in) :: zweights_SR
  complex(xp),dimension(local_nz), intent(in) :: f
  COMPLEX(xp), intent(out) :: intf
  COMPLEX(xp)              :: buffer, local_int
  INTEGER :: root, i_, iz, ierr
  ! Buil local sum using the weights of composite Simpson's rule
  local_int = 0._xp
  DO iz = 1,local_nz
      local_int = local_int + zweights_SR(iz)*f(iz)
  ENDDO
  buffer = local_int
  root = 0
  !Gather manually among the rank_z=0 processes and perform the sum
  intf = 0._xp
  IF (num_procs_z .GT. 1) THEN
      !! Everyone sends its local_sum to root = 0
      IF (rank_z .NE. root) THEN
          CALL MPI_SEND(buffer, 1 , mpi_xp_c, root, 5678, comm_z, ierr)
      ELSE
          ! Recieve from all the other processes
          DO i_ = 0,num_procs_z-1
              IF (i_ .NE. rank_z) &
                  CALL MPI_RECV(buffer, 1 , mpi_xp_c, i_, 5678, comm_z, MPI_STATUS_IGNORE, ierr)
                  intf = intf + buffer
          ENDDO
      ENDIF
      CALL manual_0D_bcast(intf)
  ELSE
    intf = local_int
  ENDIF
END SUBROUTINE simpson_rule_z

END MODULE calculus
