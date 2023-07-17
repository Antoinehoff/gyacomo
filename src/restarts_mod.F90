MODULE restarts
USE basic
USE futils,          ONLY: openf, closef, getarr, getatt, isgroup,&
                           isdataset, getarrnd, putarrnd
USE grid, ONLY: local_Na,local_Na_offset,local_np,local_np_offset,&
  local_nj,local_nj_offset,local_nky,local_nky_offset,local_nkx,local_nkx_offset,&
  local_nz,local_nz_offset,ngp,ngj,ngz, deltap,&
  total_Na,total_Np,total_Nj,total_Nky,total_Nkx,total_Nz,ieven,&
  kyarray_full,kxarray_full,zarray, zarray_full, local_zmin, local_zmax ! for z interp
USE fields
USE diagnostics_par
USE time_integration
USE prec_const, ONLY : xp,dp,sp
IMPLICIT NONE

PUBLIC :: load_moments!, write_restart

CONTAINS

  !******************************************************************************!
  !!!!!!! Fill initial moments value using the final state of a previous simulation
  !******************************************************************************!
  SUBROUTINE load_moments
    USE parallel, ONLY: comm0
    IMPLICIT NONE
    CHARACTER(LEN=50) :: dset_name
    REAL(xp):: time_cp
    INTEGER :: cstep_cp, jobnum_cp
    INTEGER :: n_
    INTEGER :: deltap_cp
    INTEGER :: n0, Np_cp, Nj_cp, Nkx_cp, Nky_cp, Nz_cp, Na_cp
    INTEGER :: ia,ip,ij,iky,ikx,iz, iacp,ipcp,ijcp,iycp,ixcp, ierr
    INTEGER :: ipi,iji,izi,izg
    REAL(xp), DIMENSION(:), ALLOCATABLE :: ky_cp, kx_cp, z_cp
    REAL(xp):: Lx_cp, Ly_cp, Lz_cp, dz_cp, zi, zj_cp, zjp1_cp
    REAL(xp):: timer_tot_1,timer_tot_2, zerotoone
    INTEGER,  DIMENSION(:), ALLOCATABLE  :: z_idx_mapping
    REAL(xp), DIMENSION(:), ALLOCATABLE  :: z_cp_stretched
    COMPLEX(xp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: moments_cp
    CALL cpu_time(timer_tot_1)
    ! Checkpoint filename
    WRITE(rstfile,'(a,a1,i2.2,a3)') TRIM(resfile0),'_',job2load,'.h5'
    CALL speak("Resume from "//rstfile)
    ! Open file
    IF(xp .EQ. dp) THEN
      CALL openf(rstfile, fidrst, mode='r', real_prec='d', mpicomm=comm0)
    ELSE
      CALL openf(rstfile, fidrst, mode='r', mpicomm=comm0)
    ENDIF
    ! Get the dimensions of the checkpoint moments
    CALL getatt(fidrst,"/data/input/model",  "Na", Na_cp)
    CALL getatt(fidrst,"/data/input/grid" ,  "Np", Np_cp)
    CALL getatt(fidrst,"/data/input/grid" ,  "Nj", Nj_cp)
    CALL getatt(fidrst,"/data/input/grid" , "Nky", Nky_cp)
    CALL getatt(fidrst,"/data/input/grid" , "Nkx", Nkx_cp)
    CALL getatt(fidrst,"/data/input/grid" ,  "Nz", Nz_cp)
    !! Check dimensional compatibility with the planned simulation
    IF(Na_cp  .NE. total_Na ) CALL speak("NOTE: Na  has changed")
    IF(Np_cp  .NE. total_Np ) CALL speak("NOTE: Np  has changed")
    IF(Nj_cp  .NE. total_Nj ) CALL speak("NOTE: Nj  has changed")
    IF(Nky_cp .NE. total_Nky) CALL speak("NOTE: Nky has changed")
    IF(Nkx_cp .NE. total_Nkx) CALL speak("NOTE: Nkx has changed")
    IF(Nz_cp  .NE. total_Nz)  CALL speak("NOTE: Nz  has changed")
    ! Get the x,y fourier modes and the z space and check grids
    ALLOCATE(ky_cp(Nky_cp),kx_cp(Nkx_cp),z_cp(Nz_cp))
    CALL getarr(fidrst,"/data/grid/coordky",ky_cp); Ly_cp = 2._xp*PI/ky_cp(2)
    CALL getarr(fidrst,"/data/grid/coordkx",kx_cp); Lx_cp = 2._xp*PI/kx_cp(2)
    CALL getarr(fidrst,"/data/grid/coordz" , z_cp); Lz_cp =-2._xp*z_cp(1)
    ! check changes
    IF(ABS(ky_cp(2)-kyarray_full(2)) .GT. 1e-3) CALL speak("NOTE: Ly  has changed")
    IF(ABS(kx_cp(2)-kxarray_full(2)) .GT. 1e-3) CALL speak("NOTE: Lx  has changed")
    IF(ABS(z_cp(1) - zarray_full(1)) .GT. 1e-3) CALL speak("NOTE: Lz  has changed")
    dz_cp = (z_cp(2)- z_cp(1)) !! check deltaz changes
    ! IF(ABS(deltaz - dz_cp) .GT. 1e-3) ERROR STOP "ERROR STOP: change of deltaz is not implemented."
    ! compute the mapping for z grid adaptation
    IF(Nz_cp .GT. 1) THEN
      ALLOCATE(z_idx_mapping(total_Nz+1),z_cp_stretched(Nz_cp+1)) ! we allocate them +1 for periodicity
      CALL z_grid_mapping(total_nz,Nz_cp,zarray_full,z_idx_mapping,z_cp_stretched)
    ENDIF
    !! 
    CALL getatt(fidrst,"/data/input/grid" ,"deltap",deltap_cp)
    IF(deltap_cp .NE. deltap) &
      ERROR STOP "!! cannot change deltap in a restart, not implemented !!"

    !!!! Looking for the latest set in the full moments output (5D arrays)
    CALL getatt(fidrst,"/data/input/basic" , "start_iframe5d", n0)
    ! Find the last results of the checkpoint file by iteration
    n_ = n0+1
    WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments", n_ ! start with moments/000001
    DO WHILE (isdataset(fidrst, dset_name)) ! If n_ is not a file we stop the loop
      n_ = n_ + 1
      WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments", n_ ! updtate file number
    ENDDO
    n_ = n_ - 1 ! n_ is not a file so take the previous one n_-1
    WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments", n_
    ! Read time dependent attributes to continue simulation
    CALL getatt(fidrst, dset_name,  'cstep', cstep_cp)
    CALL getatt(fidrst, dset_name,   'time', time_cp)
    CALL getatt(fidrst, dset_name, 'jobnum', jobnum_cp)
    ! Set the cstep, step, time and iframnd variables in basic from the checkpoint
    CALL set_basic_cp(cstep_cp,time_cp,jobnum_cp)
    CALL speak('.. restart from t = '//str(time))
    ! Read state of system from checkpoint file
    ! Brute force loading: load the full moments and take what is needed (RAM dangerous...)
    ! (other possibility would be to load all on process 0 and send the slices)
    CALL allocate_array(moments_cp, 1,Na_cp, 1,Np_cp, 1,Nj_cp, 1,Nky_cp, 1,Nkx_cp, 1,Nz_cp)
    WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments", n_
    CALL getarr(fidrst, dset_name, moments_cp)
    moments     = 0._xp;
    x: DO ikx=1,local_nkx
    ixcp = ikx+local_nkx_offset
      y: DO iky=1,local_nky
      iycp = iky + local_nky_offset
        j: DO ij=1,local_nj
        ijcp = ij + local_nj_offset
        iji  = ij + ngj/2
          p: DO ip=1,local_np
          ipcp = ip + local_np_offset
          ipi  = ip + ngp/2
            a: DO ia=1,local_na
            iacp = ia + local_na_offset
              z: DO iz = 1,local_nz
              izi     = iz + ngz/2           ! local interior index (for ghosted arrays)
              izg     = iz + local_nz_offset ! global index
              ! Checks that data exists (allows to extend the arrays if needed)
              IF((iacp.LE.Na_cp).AND.(ipcp.LE.Np_cp).AND.(ijcp.LE.Nj_cp).AND.(iycp.LE.Nky_cp).AND.(ixcp.LE.Nkx_cp)) THEN
                  IF(Nz_cp .GT. 1) THEN ! interpolate
                    zi      = zarray(izi,ieven)           ! position (=zarray_full(izg))
                    zj_cp   = z_cp_stretched(z_idx_mapping(izg))
                    zjp1_cp = z_cp_stretched(z_idx_mapping(izg)+1)
                    zerotoone = (zi - zj_cp)/(zjp1_cp-zj_cp) ! weight for interpolation
                    moments(ia,ipi,iji,iky,ikx,izi,:) = &
                      zerotoone       *moments_cp(iacp,ipcp,ijcp,iycp,ixcp,z_idx_mapping(izg))&
                    +(1._xp-zerotoone)*moments_cp(iacp,ipcp,ijcp,iycp,ixcp,z_idx_mapping(izg+1))
                  ELSE ! copy on all planes
                    moments(ia,ipi,iji,iky,ikx,izi,:) = moments_cp(iacp,ipcp,ijcp,iycp,ixcp,1)
                  ENDIF
                ELSE
                  moments(ia,ipi,iji,iky,ikx,izi,:) = 0._xp
                ENDIF
              ENDDO z
            ENDDO a
          ENDDO p
        ENDDO j
      ENDDO y
    ENDDO x
    !! deallocate the full moment variable
    DEALLOCATE(moments_cp)
    CALL closef(fidrst)
    CALL speak("Reading from restart file "//TRIM(rstfile)//" completed!")
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    ! stop time measurement
    CALL cpu_time(timer_tot_2)
    CALL speak('** Total load time : '// str(timer_tot_2 - timer_tot_1)//' **')
  END SUBROUTINE load_moments
  !******************************************************************************!
  ! Auxiliary routine
  SUBROUTINE z_grid_mapping(Nz_new,Nz_cp,zarray_new,z_idx_mapping, z_cp_stretched)
    ! Adapt check point to new zgrid when dz is kept constant
    ! We adapt the data in z by stretching and interpolating it to the new grid
    !                   o_o_o_o_o_o_o             checkpoint grid
    !                   1 2 3 4 5 6 7
    !
    ! x_____x_____x_____x_____x_____x_____x_____x new grid
    ! 1     2     3     4     5     6     7     8
    ! V                                         V (start and end are fixed)
    ! o______o______o______o______o______o______o stretched checkpoint grid (dz_s = Lz_new/Nzcp)
    ! 1      2      3      4      5      6      7 
    ! 
    ! 1_____1_____2_____3_____4_____5_____6_____7 left index mapping in z_idx_mapping
      INTEGER,  INTENT(IN) :: Nz_new, Nz_cp
      REAL(xp), DIMENSION(Nz_new),   INTENT(IN) :: zarray_new
      INTEGER,  DIMENSION(Nz_new+1), INTENT(OUT) :: z_idx_mapping
      REAL(xp), DIMENSION(Nz_cp+1) , INTENT(OUT) :: z_cp_stretched
      ! local variables
      INTEGER :: iz_new, jz_cp
      REAL(xp):: zi, zj_cp, dz_s
      LOGICAL :: in_interval
      ! stretched checkpoint grid interval 
      dz_s  = (zarray_new(Nz_new)-zarray_new(1))/(Nz_cp-1) 
      ! We loop over each new z grid points 
      DO iz_new = 1,Nz_new
        zi = zarray_new(iz_new) ! current position
        ! Loop over the stretched checkpoint grid to find the right interval
        zj_cp = zarray_new(1) ! init cp grid position
        jz_cp = 1             ! init cp grid index
        in_interval = .FALSE. ! flag to check if we stand in the interval
        DO WHILE (.NOT. in_interval)
          in_interval = (zi .GE. zj_cp) .AND. (zi .LT. zj_cp+dz_s)
          ! Increment
          zj_cp = zj_cp + dz_s
          jz_cp = jz_cp + 1
          IF(jz_cp .GT. Nz_cp+1) ERROR STOP "STOP: could not adapt grid .."
        ENDDO ! per construction the while loop should always top
        z_idx_mapping(iz_new) = jz_cp-1 ! The last index was one too much so we store the one before
      ENDDO
      ! we build explicitly the stretched cp grid for output and double check
      DO jz_cp = 1,Nz_cp
        z_cp_stretched(jz_cp) = zarray_new(1) + (jz_cp-1)*dz_s
      ENDDO
      ! Periodicity
      z_cp_stretched(Nz_cp+1) = z_cp_stretched(1)
      z_idx_mapping (Nz_new+1) = z_idx_mapping (1)
      ! Check that the start and the end of the grids are the same
      IF(.NOT.(abs(z_cp_stretched(1)-zarray_new(1)) .LT. 1e-3).AND.(abs(z_cp_stretched(Nz_cp)-zarray_new(Nz_new)).LT.1e-3)) &
        ERROR STOP "Failed to stretch the cp grid"
    END SUBROUTINE z_grid_mapping
END MODULE restarts
