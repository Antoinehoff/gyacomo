MODULE utility

  USE basic

  use prec_const
  IMPLICIT NONE
  PUBLIC :: manual_2D_bcast, manual_3D_bcast,&
            simpson_rule_z, o2e_z, e2o_z

CONTAINS

  FUNCTION is_nan(x,str) RESULT(isn)

    USE time_integration

    use prec_const
    IMPLICIT NONE

    real(dp), INTENT(IN) :: x
    CHARACTER(LEN=*), INTENT(IN) :: str
    LOGICAL :: isn

    isn=.FALSE.
    IF(x .NE. x) THEN
       isn=.TRUE.
    END IF
    IF((isn).AND.(str.NE.'')) THEN
       WRITE(*,'(a20,a25,i6.6,a20,i1)') str,' = NaN  at timestep',cstep, ' and substep',updatetlevel
       CALL FLUSH(stdout)
    END IF
  END FUNCTION is_nan


  FUNCTION is_inf(x,str) RESULT(isi)

    USE time_integration

    use prec_const
    IMPLICIT NONE

    real(dp), INTENT(IN) :: x
    CHARACTER(LEN=*), INTENT(IN) :: str
    LOGICAL :: isi

    isi=.FALSE.
    IF(x+1.0== x) THEN
       isi=.TRUE.
    END IF
    IF((isi).AND.(str.NE.'')) THEN
       !WRITE(*,'(a20,a25,i6.6,a20,i1)') str,' = Inf at timestep',cstep, ' and substep',updatetlevel
       CALL FLUSH(stdout)
    END IF

  END FUNCTION is_inf

  FUNCTION checkfield(field,str) RESULT(mlend)

    USE grid

    use prec_const
    IMPLICIT NONE

    COMPLEX(dp), DIMENSION(ikxs:ikxe,ikys:ikye), INTENT(IN) :: field
    CHARACTER(LEN=*), INTENT(IN) :: str
    LOGICAL :: mlend
    COMPLEX(dp) :: sumfield

    sumfield=SUM(field)

    mlend= is_nan( REAL(sumfield),str).OR.is_inf( REAL(sumfield),str) &
      .OR. is_nan(AIMAG(sumfield),str).OR.is_inf(AIMAG(sumfield),str)
  END FUNCTION checkfield


    !!!!! This is a manual way to do MPI_BCAST !!!!!!!!!!!
  SUBROUTINE manual_2D_bcast(field_)
    USE grid
    IMPLICIT NONE
    COMPLEX(dp), INTENT(INOUT) :: field_(ikxs:ikxe,ikys:ikye)
    COMPLEX(dp) :: buffer(ikxs:ikxe,ikys:ikye)
    INTEGER     :: i_, root, world_rank, world_size
    root = 0;
    CALL MPI_COMM_RANK(comm_p,world_rank,ierr)
    CALL MPI_COMM_SIZE(comm_p,world_size,ierr)
    IF (world_size .GT. 1) THEN
      !! Broadcast phi to the other processes on the same k range (communicator along p)
      IF (world_rank .EQ. root) THEN
        ! Fill the buffer
        DO ikx = ikxs,ikxe
          DO iky = ikys,ikye
            buffer(ikx,iky) = field_(ikx,iky)
          ENDDO
        ENDDO
        ! Send it to all the other processes
        DO i_ = 0,num_procs_p-1
          IF (i_ .NE. world_rank) &
          CALL MPI_SEND(buffer, local_nkx * Nky , MPI_DOUBLE_COMPLEX, i_, 0, comm_p, ierr)
        ENDDO
      ELSE
        ! Recieve buffer from root
        CALL MPI_RECV(buffer, local_nkx * Nky , MPI_DOUBLE_COMPLEX, root, 0, comm_p, MPI_STATUS_IGNORE, ierr)
        ! Write it in phi
        DO ikx = ikxs,ikxe
          DO iky = ikys,ikye
            field_(ikx,iky) = buffer(ikx,iky)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  END SUBROUTINE manual_2D_bcast

  !!!!! This is a manual way to do MPI_BCAST !!!!!!!!!!!
SUBROUTINE manual_3D_bcast(field_)
  USE grid
  IMPLICIT NONE
  COMPLEX(dp), INTENT(INOUT) :: field_(ikxs:ikxe,ikys:ikye,izs:ize)
  COMPLEX(dp) :: buffer(ikxs:ikxe,ikys:ikye,izs:ize)
  INTEGER     :: i_, root, world_rank, world_size
  root = 0;
  CALL MPI_COMM_RANK(comm_p,world_rank,ierr)
  CALL MPI_COMM_SIZE(comm_p,world_size,ierr)
  IF (world_size .GT. 1) THEN
    !! Broadcast phi to the other processes on the same k range (communicator along p)
    IF (world_rank .EQ. root) THEN
      ! Fill the buffer
      DO ikx = ikxs,ikxe
        DO iky = ikys,ikye
          DO iz = izs,ize
            buffer(ikx,iky,iz) = field_(ikx,iky,iz)
          ENDDO
        ENDDO
      ENDDO
      ! Send it to all the other processes
      DO i_ = 0,num_procs_p-1
        IF (i_ .NE. world_rank) &
        CALL MPI_SEND(buffer, local_nkx * Nky * Nz, MPI_DOUBLE_COMPLEX, i_, 0, comm_p, ierr)
      ENDDO
    ELSE
      ! Recieve buffer from root
      CALL MPI_RECV(buffer, local_nkx * Nky * Nz, MPI_DOUBLE_COMPLEX, root, 0, comm_p, MPI_STATUS_IGNORE, ierr)
      ! Write it in phi
      DO ikx = ikxs,ikxe
        DO iky = ikys,ikye
          DO iz = izs,ize
            field_(ikx,iky,iz) = buffer(ikx,iky,iz)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF
END SUBROUTINE manual_3D_bcast

SUBROUTINE simpson_rule_z(f,intf)
 ! integrate f(z) over z using the simpon's rule. Assume periodic boundary conditions (f(ize+1) = f(izs))
 !from molix BJ Frei
 use prec_const
 use grid
 !
 implicit none
 !
 complex(dp),dimension(izs:ize), intent(in) :: f
 COMPLEX(dp), intent(out) :: intf
 !
 COMPLEX(dp) :: buff_
 !
 IF(Nz .GT. 1) THEN
   IF(mod(Nz,2) .ne. 0 ) THEN
      ERROR STOP 'Simpson rule: Nz must be an even number  !!!!'
   ENDIF
   !
   buff_ = 0._dp
   !
   DO iz = izs, Nz/2
      IF(iz .eq. Nz/2) THEN ! ... iz = ize
         buff_ = buff_ + (f(izs) + 4._dp*f(ize) + f(ize-1 ))
      ELSE
         buff_ = buff_ + (f(2*iz+1) + 4._dp*f(2*iz) + f(2*iz-1 ))
      ENDIF
   ENDDO
   !
   !
   intf = buff_*deltaz/3._dp
   !
 ELSE
   intf = f(izs)
 ENDIF
END SUBROUTINE simpson_rule_z

SUBROUTINE o2e_z(fo,fe)
 ! gives the value of a field from the odd grid to the even one
 use prec_const
 use grid
 !
 implicit none
 !
 COMPLEX(dp),dimension(1:Nz), intent(in)  :: fo
 COMPLEX(dp),dimension(1:Nz), intent(out) :: fe !
 !
 DO iz = 2,Nz
   fe(iz) = 0.5_dp*(fo(iz)+fo(iz-1))
 ENDDO
 ! Periodic boundary conditions
 fe(1) = 0.5_dp*(fo(1) + fo(Nz))
END SUBROUTINE o2e_z

SUBROUTINE e2o_z(fe,fo)
 ! gives the value of a field from the even grid to the odd one
 use prec_const
 use grid
 !
 implicit none
 !
 COMPLEX(dp),dimension(1:Nz), intent(in)  :: fe
 COMPLEX(dp),dimension(1:Nz), intent(out) :: fo
 !
 DO iz = 1,Nz-1
   fo(iz) = 0.5_dp*(fe(iz+1)+fe(iz))
 ENDDO
 ! Periodic boundary conditions
 fo(Nz) = 0.5_dp*(fe(1) + fe(Nz))
END SUBROUTINE e2o_z

END MODULE utility
