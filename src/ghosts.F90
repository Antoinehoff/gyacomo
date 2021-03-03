
!Communicate two extra buffer zones needed for using left/right fd schemes in parallel 
!gradients routines
SUBROUTINE update_ghosts_p
    USE basic
    USE fields
    USE grid
    USE ppinit
  
    use prec_const
    IMPLICIT NONE
     
    INTEGER :: status(MPI_STATUS_SIZE)
    complex(dp):: buff_snd_L( 1: 2,ijs_e:ije_e,ikrs:ikre,ikzs:ikze)
    complex(dp):: buff_snd_R(-2:-1,ijs_e:ije_e,ikrs:ikre,ikzs:ikze)
    complex(dp):: buff_rcv_R( 1: 2,ijs_e:ije_e,ikrs:ikre,ikzs:ikze)
    complex(dp):: buff_rcv_L(-2:-1,ijs_e:ije_e,ikrs:ikre,ikzs:ikze)
  
    !! Set up data to send to left neighbor
    DO ij=ijs_e,ije_e
        DO ikr=ikrs,ikre
            DO ikz=ikzs,ikze
                ! send to left ipe + 1 moments
                buff_snd_L(1,ij,ikr,ikz) = moments_e(ips,ij,ikr,ikz,updatetlevel)
                ! send to left ipe + 2 moments
                buff_snd_L(2,ij,ikr,ikz) = moments_e(ips+1,ij,ikr,ikz,updatetlevel)
            END DO
        END DO
    END DO

    ! Exchange data with left neighbor
    CALL mpi_sendrecv(buff_snd_L, count, MPI_DOUBLE_PRECISION, left_neighbor, 0, &
    buff_rcv_L, count, MPI_DOUBLE_PRECISION, source, 0, &
    commx, status, ierr)
   
    ! Write received data
    DO ij=ijs_e,ije_e
        DO ikr=ikrs,ikre
            DO ikz=ikzs,ikze  
                ! ip - 1 moments
                moments_e(ips-1,ij,ikr,ikz,updatetlevel) = buff_rcv_L(-1,ij,ikr,ikz)
                ! ip - 2 moments
                moments_e(ips-2,ij,ikr,ikz,updatetlevel) = buff_rcv_L(-2,ij,ikr,ikz)
            END DO
        END DO
    END DO 

    !! Set up data to send to right neighbors
    DO ij=ijs_e,ije_e
        DO ikr=ikrs,ikre
            DO ikz=ikzs,ikze
                ! send to right ips - 1 moments
                buff_snd_R(-1,ij,ikr,ikz) = moments_e(ipe,ij,ikr,ikz,updatetlevel)
                ! send to right ips - 2 moments
                buff_snd_R(-2,ij,ikr,ikz) = moments_e(ipe-1,ij,ikr,ikz,updatetlevel)
            END DO
        END DO
    END DO
   
    ! Exchange data with right neighbor
    CALL mpi_sendrecv(buff_snd_R, count, MPI_DOUBLE_PRECISION, dest, 0, &
    buff_rcv_R, count, MPI_DOUBLE_PRECISION, source, 0, &
    commx, status, ierr)

    ! Write received data
    DO ij=ijs_e,ije_e
        DO ikr=ikrs,ikre
            DO ikz=ikzs,ikze  
                ! ipe + 1 moments
                moments_e(ipe+1,ij,ikr,ikz,updatetlevel) = buff_rcv_R(1,ij,ikr,ikz)
                ! ip + 2 moments
                moments_e(ipe+2,ij,ikr,ikz,updatetlevel) = buff_rcv_R(2,ij,ikr,ikz)
            END DO
        END DO
    END DO 

    
  END SUBROUTINE update_ghosts_p
  