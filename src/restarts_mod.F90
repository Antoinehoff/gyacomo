MODULE restarts
USE basic
USE futils,          ONLY: openf, closef, getarr, getatt, isgroup,&
                           isdataset, getarrnd, putarrnd
USE grid
USE fields
USE diagnostics_par
USE time_integration
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
    INTEGER :: pmax_cp, jmax_cp, n0, Nkx_cp, Nky_cp, Nz_cp, Na_cp, Np_cp, Nj_cp
    INTEGER :: ia,ip,ij,iky,ikx,iz, iacp,ipcp,ijcp,iycp,ixcp,izcp, ierr
    INTEGER :: ipi,iji,izi
    REAL(xp):: timer_tot_1,timer_tot_2
    COMPLEX(xp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: moments_cp
    CALL cpu_time(timer_tot_1)
    ! Checkpoint filename
    WRITE(rstfile,'(a,a1,i2.2,a3)') TRIM(resfile0),'_',job2load,'.h5'
    CALL speak("Resume from "//rstfile)
    ! Open file
    CALL openf(rstfile, fidrst,mpicomm=comm0)
    ! Get the checkpoint moments degrees to allocate memory
    CALL getatt(fidrst,"/data/input/grid" ,   "Nkx",   Nkx_cp)
    CALL getatt(fidrst,"/data/input/grid" ,   "Nky",   Nky_cp)
    CALL getatt(fidrst,"/data/input/grid" ,    "Nz",    Nz_cp)
    IF(Nz_cp .NE. Nz) &
      ERROR STOP "!! cannot change Nz in a restart, interp or reduction not implemented !!"
    CALL getatt(fidrst,"/data/input/grid" ,"deltap",deltap_cp)
    IF(deltap_cp .NE. deltap) &
      ERROR STOP "!! cannot change deltap in a restart, not implemented !!"
    CALL getatt(fidrst,"/data/input/grid" , "pmax", pmax_cp)
    Np_cp = pmax_cp/deltap_cp+1
    CALL getatt(fidrst,"/data/input/grid" , "jmax", jmax_cp)
    Nj_cp = jmax_cp+1
    CALL getatt(fidrst,"/data/input/model",   "Na",   Na_cp)
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
    ! other possibility is to loop over slices
    CALL allocate_array(moments_cp, 1,Na_cp, 1,Np_cp, 1,Nj_cp, 1,Nky_cp, 1,Nkx_cp, 1,Nz_cp)
    WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments", n_
    CALL getarr(fidrst, dset_name, moments_cp(:,:,:,:,:,:))

    moments     = 0._xp;
    z: DO iz = 1,local_nz
      izcp = iz + local_nz_offset
      izi  = iz + ngz/2
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
              a: DO ia=1,Na_cp
                iacp = ia + local_na_offset
                ! IF((iacp.LE.Na_cp).AND.(ipcp.LE.Np_cp).AND.(ijcp.LE.Nj_cp).AND.(iycp.LE.Nky_cp).AND.(ixcp.LE.Nkx_cp).AND.(izcp.LE.Nz_cp)) &
                  moments(ia,ipi,iji,iky,ikx,izi,1) = moments_cp(iacp,ipcp,ijcp,iycp,ixcp,izcp)
              ENDDO a
            ENDDO p
          ENDDO j
        ENDDO y
      ENDDO x
    ENDDO z
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

END MODULE restarts
