SUBROUTINE diagnose(kstep)
  !   Diagnostics, writing simulation state to disk

  USE basic
  USE grid
  USE diagnostics_par
  USE futils, ONLY: creatf, creatg, creatd, closef, putarr, putfile, attach, openf, putarrnd
  USE model
  USE initial_par
  USE fields
  USE time_integration
  USE utility
  USE prec_const
  IMPLICIT NONE

  INCLUDE 'srcinfo.h'

  INTEGER, INTENT(in) :: kstep
  INTEGER, parameter  :: BUFSIZE = 2
  INTEGER :: rank = 0
  INTEGER :: dims(1) = (/0/)
  INTEGER :: cp_counter = 0
  CHARACTER(len=256) :: str, fname,test_

  CALL cpu_time(t0_diag) ! Measuring time
  !_____________________________________________________________________________
  !                   1.   Initial diagnostics

  IF ((kstep .EQ. 0)) THEN
    ! Writing output filename
     WRITE(resfile,'(a,a1,i2.2,a3)') TRIM(resfile0),'_',jobnum,'.h5'
     !                      1.1   Initial run
     ! Main output file creation
     IF (write_doubleprecision) THEN
        CALL creatf(resfile, fidres, real_prec='d', mpicomm=comm0)
     ELSE
        CALL creatf(resfile, fidres, mpicomm=comm0)
     END IF
     IF (my_id .EQ. 0) WRITE(*,'(3x,a,a)') TRIM(resfile), ' created'

     !  Data group
     CALL creatg(fidres, "/data", "data")

     !  File group
     CALL creatg(fidres, "/files", "files")
     CALL attach(fidres, "/files",  "jobnum", jobnum)

     ! Profiler time measurement
     CALL creatg(fidres, "/profiler", "performance analysis")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_rhs",        "cumulative rhs computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_adv_field",  "cumulative adv. fields computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_clos",       "cumulative closure computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_ghost",       "cumulative communication time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_coll",       "cumulative collision computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_poisson",    "cumulative poisson computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_Sapj",       "cumulative Sapj computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_checkfield", "cumulative checkfield computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_diag",       "cumulative sym computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_step",       "cumulative total step computation time")
     CALL creatd(fidres, 0, dims, "/profiler/time",          "current simulation time")

     ! Grid info
     CALL creatg(fidres, "/data/grid", "Grid data")
     CALL putarr(fidres, "/data/grid/coordkx",   kxarray_full,  "kx*rho_s0", ionode=0)
     CALL putarr(fidres, "/data/grid/coordky",   kyarray_full,  "ky*rho_s0", ionode=0)
     CALL putarr(fidres, "/data/grid/coordz",    zarray_full,   "z/R", ionode=0)
     CALL putarr(fidres, "/data/grid/coordp_e" , parray_e_full, "p_e", ionode=0)
     CALL putarr(fidres, "/data/grid/coordj_e" , jarray_e_full, "j_e", ionode=0)
     CALL putarr(fidres, "/data/grid/coordp_i" , parray_i_full, "p_i", ionode=0)
     CALL putarr(fidres, "/data/grid/coordj_i" , jarray_i_full, "j_i", ionode=0)

     !  var0d group (gyro transport)
     IF (nsave_0d .GT. 0) THEN
      CALL creatg(fidres, "/data/var0d", "0d profiles")
      CALL creatd(fidres, rank, dims,  "/data/var0d/time",     "Time t*c_s/R")
      CALL creatd(fidres, rank, dims, "/data/var0d/cstep", "iteration number")

      IF (write_gamma) THEN
        CALL creatd(fidres, rank, dims, "/data/var0d/gflux_ri", "Radial gyro ion transport")
        CALL creatd(fidres, rank, dims, "/data/var0d/pflux_ri", "Radial part ion transport")
      ENDIF
      IF (write_hf) THEN
        CALL creatd(fidres, rank, dims, "/data/var0d/hflux_x", "Radial part ion heat flux")
      ENDIF
      IF (cstep==0) THEN
        iframe0d=0
      ENDIF
      CALL attach(fidres,"/data/var0d/" , "frames", iframe0d)
     END IF


     !  var2d group (??)
     IF (nsave_2d .GT. 0) THEN
      CALL creatg(fidres, "/data/var2d", "2d profiles")
      CALL creatd(fidres, rank, dims,  "/data/var2d/time",     "Time t*c_s/R")
      CALL creatd(fidres, rank, dims, "/data/var2d/cstep", "iteration number")
      IF (cstep==0) THEN
        iframe2d=0
      ENDIF
      CALL attach(fidres,"/data/var2d/" , "frames", iframe2d)
     END IF

     !  var3d group (electro. pot., Ni00 moment)
     IF (nsave_3d .GT. 0) THEN
      CALL creatg(fidres, "/data/var3d", "3d profiles")
      CALL creatd(fidres, rank, dims,  "/data/var3d/time",     "Time t*c_s/R")
      CALL creatd(fidres, rank, dims, "/data/var3d/cstep", "iteration number")

      IF (write_phi) CALL creatg(fidres, "/data/var3d/phi", "phi")

      IF (write_Na00) THEN
       CALL creatg(fidres, "/data/var3d/Ne00", "Ne00")
       CALL creatg(fidres, "/data/var3d/Ni00", "Ni00")
      ENDIF

      IF (write_dens) THEN
       CALL creatg(fidres, "/data/var3d/dens_e", "dens_e")
       CALL creatg(fidres, "/data/var3d/dens_i", "dens_i")
      ENDIF

      IF (write_temp) THEN
       CALL creatg(fidres, "/data/var3d/temp_e", "temp_e")
       CALL creatg(fidres, "/data/var3d/temp_i", "temp_i")
      ENDIF

      IF (cstep==0) THEN
        iframe3d=0
      ENDIF
      CALL attach(fidres,"/data/var3d/" , "frames", iframe3d)
     END IF

     !  var5d group (moments)
     IF (nsave_5d .GT. 0) THEN
       CALL creatg(fidres, "/data/var5d", "5d profiles")
       CALL creatd(fidres, rank, dims,  "/data/var5d/time",     "Time t*c_s/R")
       CALL creatd(fidres, rank, dims, "/data/var5d/cstep", "iteration number")

       IF (write_Napj) THEN
        CALL creatg(fidres, "/data/var5d/moments_e", "moments_e")
        CALL creatg(fidres, "/data/var5d/moments_i", "moments_i")
       ENDIF

       IF (write_Sapj) THEN
        CALL creatg(fidres, "/data/var5d/moments_e", "Sipj")
        CALL creatg(fidres, "/data/var5d/moments_i", "Sepj")
       ENDIF

       IF (cstep==0) THEN
        iframe5d=0
       END IF
       CALL attach(fidres,"/data/var5d/" , "frames", iframe5d)
     END IF

     !  Add input namelist variables as attributes of /data/input, defined in srcinfo.h
     IF (my_id .EQ. 0) WRITE(*,*) 'VERSION=', VERSION
     IF (my_id .EQ. 0) WRITE(*,*)  'BRANCH=', BRANCH
     IF (my_id .EQ. 0) WRITE(*,*)  'AUTHOR=', AUTHOR
     IF (my_id .EQ. 0) WRITE(*,*)    'HOST=', HOST

     WRITE(str,'(a,i2.2)') "/data/input"
     CALL creatd(fidres, rank,dims,TRIM(str),'Input parameters')
     CALL attach(fidres, TRIM(str),     "version",  VERSION) !defined in srcinfo.h
     CALL attach(fidres, TRIM(str),      "branch",   BRANCH) !defined in srcinfo.h
     CALL attach(fidres, TRIM(str),      "author",   AUTHOR) !defined in srcinfo.h
     CALL attach(fidres, TRIM(str),    "execdate", EXECDATE) !defined in srcinfo.h
     CALL attach(fidres, TRIM(str),        "host",     HOST) !defined in srcinfo.h
     CALL attach(fidres, TRIM(str),  "start_time",     time)
     CALL attach(fidres, TRIM(str), "start_cstep",    cstep-1)
     CALL attach(fidres, TRIM(str), "start_iframe0d", iframe0d)
     CALL attach(fidres, TRIM(str), "start_iframe2d", iframe2d)
     CALL attach(fidres, TRIM(str), "start_iframe3d", iframe3d)
     CALL attach(fidres, TRIM(str), "start_iframe5d", iframe5d)
     CALL attach(fidres, TRIM(str),          "dt",       dt)
     CALL attach(fidres, TRIM(str),        "tmax",     tmax)
     CALL attach(fidres, TRIM(str),        "nrun",     nrun)
     CALL attach(fidres, TRIM(str),    "cpu_time",       -1)
     CALL attach(fidres, TRIM(str),       "Nproc",   num_procs)
     CALL attach(fidres, TRIM(str),       "Np_p" , num_procs_p)
     CALL attach(fidres, TRIM(str),       "Np_kx",num_procs_kx)
     CALL attach(fidres, TRIM(str), "write_gamma", write_gamma)
     CALL attach(fidres, TRIM(str),    "write_hf",    write_hf)
     CALL attach(fidres, TRIM(str),   "write_phi",   write_phi)
     CALL attach(fidres, TRIM(str),  "write_Na00",  write_Na00)
     CALL attach(fidres, TRIM(str),  "write_Napj",  write_Napj)
     CALL attach(fidres, TRIM(str),  "write_Sapj",  write_Sapj)
     CALL attach(fidres, TRIM(str),  "write_dens",  write_dens)
     CALL attach(fidres, TRIM(str),  "write_temp",  write_temp)

     CALL grid_outputinputs(fidres, str)

     CALL diag_par_outputinputs(fidres, str)

     CALL model_outputinputs(fidres, str)

     CALL initial_outputinputs(fidres, str)

     CALL time_integration_outputinputs(fidres, str)


     !  Save STDIN (input file) of this run
     IF(jobnum .LE. 99) THEN
        WRITE(str,'(a,i2.2)') "/files/STDIN.",jobnum
     ELSE
        WRITE(str,'(a,i3.2)') "/files/STDIN.",jobnum
     END IF

     INQUIRE(unit=lu_in, name=fname)
     CLOSE(lu_in)

     CALL putfile(fidres, TRIM(str), TRIM(fname),ionode=0)

   ELSEIF((kstep .EQ. 0)) THEN

      IF(jobnum .LE. 99) THEN
         WRITE(resfile,'(a,a1,i2.2,a3)') TRIM(resfile0),'_',jobnum,'.h5'
      ELSE
         WRITE(resfile,'(a,a1,i3.2,a3)') TRIM(resfile0),'_',jobnum,'.h5'
      END IF

      CALL openf(resfile,fidres, 'D');

   ENDIF

  !_____________________________________________________________________________
  !                   2.   Periodic diagnostics
  !
  IF (kstep .GE. 0) THEN

    ! Terminal info
    IF (MOD(cstep, INT(1.0/dt)) == 0 .AND. (my_id .EQ. 0)) THEN
     WRITE(*,"(F6.0,A,F6.0)") time,"/",tmax
    ENDIF

     !                       2.1   0d history arrays
     IF (nsave_0d .GT. 0) THEN
        IF ( MOD(cstep, nsave_0d) == 0 ) THEN
           CALL diagnose_0d
        END IF
     END IF

     !                       2.2   1d profiles
     ! empty in our case

     !                       2.3   2d profiles
     ! empty in our case


     !                       2.3   2d profiles
     IF (nsave_3d .GT. 0) THEN
        IF (MOD(cstep, nsave_3d) == 0) THEN
           CALL diagnose_3d
        END IF
     END IF

     !                       2.4   3d profiles
     IF (nsave_5d .GT. 0 .AND. cstep .GT. 0) THEN
        IF (MOD(cstep, nsave_5d) == 0) THEN
           CALL diagnose_5d
        END IF
     END IF

  !_____________________________________________________________________________
  !                   3.   Final diagnostics

  ELSEIF (kstep .EQ. -1) THEN
     CALL cpu_time(finish)
     CALL attach(fidres, "/data/input","cpu_time",finish-start)

     ! make a checkpoint at last timestep if not crashed
     IF(.NOT. crashed) THEN
       IF(my_id .EQ. 0) write(*,*) 'Saving last state'
       CALL diagnose_5d
     ENDIF
     ! Display computational time cost
     IF (my_id .EQ. 0) CALL display_h_min_s(finish-start)

     !   Close all diagnostic files
     CALL closef(fidres)

  END IF

  CALL cpu_time(t1_diag); tc_diag = tc_diag + (t1_diag - t0_diag)

END SUBROUTINE diagnose


SUBROUTINE diagnose_0d

  USE basic
  USE futils, ONLY: append, attach, getatt
  USE diagnostics_par
  USE prec_const
  USE processing

  IMPLICIT NONE
  ! Time measurement data
  CALL append(fidres, "/profiler/Tc_rhs",              tc_rhs,ionode=0)
  CALL append(fidres, "/profiler/Tc_adv_field",  tc_adv_field,ionode=0)
  CALL append(fidres, "/profiler/Tc_clos",            tc_clos,ionode=0)
  CALL append(fidres, "/profiler/Tc_ghost",          tc_ghost,ionode=0)
  CALL append(fidres, "/profiler/Tc_coll",            tc_coll,ionode=0)
  CALL append(fidres, "/profiler/Tc_poisson",      tc_poisson,ionode=0)
  CALL append(fidres, "/profiler/Tc_Sapj",            tc_Sapj,ionode=0)
  CALL append(fidres, "/profiler/Tc_checkfield",tc_checkfield,ionode=0)
  CALL append(fidres, "/profiler/Tc_diag",            tc_diag,ionode=0)
  CALL append(fidres, "/profiler/Tc_step",            tc_step,ionode=0)
  CALL append(fidres, "/profiler/time",                  time,ionode=0)
  ! Processing data
  CALL append(fidres,  "/data/var0d/time",           time,ionode=0)
  CALL append(fidres, "/data/var0d/cstep", real(cstep,dp),ionode=0)
  CALL getatt(fidres,      "/data/var0d/",       "frames",iframe2d)
  iframe0d=iframe0d+1
  CALL attach(fidres,"/data/var0d/" , "frames", iframe0d)
  ! Ion transport data
  IF (write_gamma) THEN
    CALL compute_radial_ion_transport
    CALL append(fidres, "/data/var0d/gflux_ri",gflux_ri,ionode=0)
    CALL append(fidres, "/data/var0d/pflux_ri",pflux_ri,ionode=0)
  ENDIF
  IF (write_hf) THEN
    CALL compute_radial_heatflux
    CALL append(fidres, "/data/var0d/hflux_x",hflux_x,ionode=0)
  ENDIF
END SUBROUTINE diagnose_0d


SUBROUTINE diagnose_2d

  USE basic
  USE futils, ONLY: append, getatt, attach, putarrnd
  USE fields
  USE array, ONLY: Ne00, Ni00, dens_e, dens_i, temp_e, temp_i
  USE grid, ONLY: ikxs,ikxe, ikys,ikye, Nkx, Nky, local_nkx, ikx, iky, ips_e, ips_i
  USE time_integration
  USE diagnostics_par
  USE prec_const
  USE processing

  IMPLICIT NONE

  COMPLEX(dp) :: buffer(ikxs:ikxe,ikys:ikye)
  INTEGER     :: i_, root, world_rank, world_size

  CALL append(fidres,  "/data/var2d/time",           time,ionode=0)
  CALL append(fidres, "/data/var2d/cstep", real(cstep,dp),ionode=0)
  CALL getatt(fidres,      "/data/var2d/",       "frames",iframe2d)
  iframe2d=iframe2d+1
  CALL attach(fidres,"/data/var2d/" , "frames", iframe2d)

CONTAINS

  SUBROUTINE write_field2d(field, text)
    USE futils, ONLY: attach, putarr
    USE grid, ONLY: ikxs,ikxe, ikys,ikye, Nkx, Nky, local_nkx
    USE prec_const
    USE basic, ONLY : comm_kx, num_procs_p, rank_p
    IMPLICIT NONE

    COMPLEX(dp), DIMENSION(ikxs:ikxe, ikys:ikye), INTENT(IN) :: field
    CHARACTER(*), INTENT(IN) :: text
    COMPLEX(dp) :: buffer_dist(ikxs:ikxe,ikys:ikye)
    COMPLEX(dp) :: buffer_full(1:Nkx,1:Nky)
    INTEGER     :: scount, rcount
    CHARACTER(LEN=50) :: dset_name

    scount = (ikxe-ikxs+1) * (ikye-ikys+1)
    rcount = scount

    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var2d", TRIM(text), iframe2d
    IF (num_procs .EQ. 1) THEN ! no data distribution
      CALL putarr(fidres, dset_name, field(ikxs:ikxe, ikys:ikye), ionode=0)
    ELSE
      CALL putarrnd(fidres, dset_name, field(ikxs:ikxe, ikys:ikye),  (/1, 1/))
    ENDIF
    CALL attach(fidres, dset_name, "time", time)

  END SUBROUTINE write_field2d

END SUBROUTINE diagnose_2d

SUBROUTINE diagnose_3d

  USE basic
  USE futils, ONLY: append, getatt, attach, putarrnd
  USE fields
  USE array, ONLY: Ne00, Ni00, dens_e, dens_i, temp_e, temp_i
  USE grid, ONLY: ikxs,ikxe, ikys,ikye, Nkx, Nky, local_nkx, ikx, iky, ips_e, ips_i
  USE time_integration
  USE diagnostics_par
  USE prec_const
  USE processing

  IMPLICIT NONE

  INTEGER     :: i_, root, world_rank, world_size

  CALL append(fidres,  "/data/var3d/time",           time,ionode=0)
  CALL append(fidres, "/data/var3d/cstep", real(cstep,dp),ionode=0)
  CALL getatt(fidres,      "/data/var3d/",       "frames",iframe3d)
  iframe3d=iframe3d+1
  CALL attach(fidres,"/data/var3d/" , "frames", iframe3d)

  IF (write_phi) CALL write_field3d(phi (:,:,:), 'phi')

  IF (write_Na00) THEN
    IF ( (ips_e .EQ. 1) .AND. (ips_i .EQ. 1) ) THEN
      Ne00(ikxs:ikxe,ikys:ikye,izs:ize) = moments_e(ips_e,1,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel)
      Ni00(ikxs:ikxe,ikys:ikye,izs:ize) = moments_i(ips_e,1,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel)
    ENDIF

    CALL manual_3D_bcast(Ne00(ikxs:ikxe,ikys:ikye,izs:ize))
    CALL write_field3d(Ne00(ikxs:ikxe,ikys:ikye,izs:ize), 'Ne00')

    CALL manual_3D_bcast(Ni00(ikxs:ikxe,ikys:ikye,izs:ize))
    CALL write_field3d(Ni00(ikxs:ikxe,ikys:ikye,izs:ize), 'Ni00')
  ENDIF

  IF (write_dens) THEN
    CALL compute_density
    CALL write_field3d(dens_e(ikxs:ikxe,ikys:ikye,izs:ize), 'dens_e')
    CALL write_field3d(dens_i(ikxs:ikxe,ikys:ikye,izs:ize), 'dens_i')
  ENDIF

  IF (write_temp) THEN
    CALL compute_temperature
    CALL write_field3d(temp_e(ikxs:ikxe,ikys:ikye,izs:ize), 'temp_e')
    CALL write_field3d(temp_i(ikxs:ikxe,ikys:ikye,izs:ize), 'temp_i')
  ENDIF

CONTAINS

  SUBROUTINE write_field3d(field, text)
    USE futils, ONLY: attach, putarr
    USE grid, ONLY: ikxs,ikxe, ikys,ikye, Nkx, Nky, local_nkx
    USE prec_const
    USE basic, ONLY : comm_kx, num_procs_p, rank_p
    IMPLICIT NONE

    COMPLEX(dp), DIMENSION(ikxs:ikxe, ikys:ikye, izs:ize), INTENT(IN) :: field
    CHARACTER(*), INTENT(IN) :: text
    CHARACTER(LEN=50) :: dset_name

    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var3d", TRIM(text), iframe3d
    IF (num_procs .EQ. 1) THEN ! no data distribution
      CALL putarr(fidres, dset_name, field(ikxs:ikxe, ikys:ikye, izs:ize), ionode=0)
    ELSE
      CALL putarrnd(fidres, dset_name, field(ikxs:ikxe, ikys:ikye, izs:ize),  (/1, 1/))
    ENDIF
    CALL attach(fidres, dset_name, "time", time)

  END SUBROUTINE write_field3d

END SUBROUTINE diagnose_3d

SUBROUTINE diagnose_5d

   USE basic
   USE futils, ONLY: append, getatt, attach, putarrnd
   USE fields
   USE array!, ONLY: Sepj, Sipj
   USE grid, ONLY: ips_e,ipe_e, ips_i, ipe_i, &
                   ijs_e,ije_e, ijs_i, ije_i
   USE time_integration
   USE diagnostics_par
   USE prec_const
   IMPLICIT NONE

   CALL append(fidres,  "/data/var5d/time",           time,ionode=0)
   CALL append(fidres, "/data/var5d/cstep", real(cstep,dp),ionode=0)
   CALL getatt(fidres,      "/data/var5d/",       "frames",iframe5d)
   iframe5d=iframe5d+1
   CALL attach(fidres,"/data/var5d/" , "frames", iframe5d)

   IF (write_Napj) THEN
    CALL write_field5d_e(moments_e(ips_e:ipe_e,ijs_e:ije_e,:,:,:,updatetlevel), 'moments_e')
    CALL write_field5d_i(moments_i(ips_i:ipe_i,ijs_i:ije_i,:,:,:,updatetlevel), 'moments_i')
   ENDIF

   IF (write_Sapj) THEN
     CALL write_field5d_e(Sepj(ips_e:ipe_e,ijs_e:ije_e,:,:,:), 'Sepj')
     CALL write_field5d_i(Sipj(ips_i:ipe_i,ijs_i:ije_i,:,:,:), 'Sipj')
   ENDIF

 CONTAINS

   SUBROUTINE write_field5d_e(field, text)
     USE futils, ONLY: attach, putarr, putarrnd
     USE grid,   ONLY: ips_e,ipe_e, ijs_e,ije_e, ikxs,ikxe, ikys,ikye, izs,ize
     USE prec_const
     IMPLICIT NONE

     COMPLEX(dp), DIMENSION(ips_e:ipe_e,ijs_e:ije_e,ikxs:ikxe,ikys:ikye,izs:ize), INTENT(IN) :: field
     CHARACTER(*), INTENT(IN) :: text

     CHARACTER(LEN=50) :: dset_name

     WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var5d", TRIM(text), iframe5d
     IF (num_procs .EQ. 1) THEN
       CALL putarr(fidres, dset_name, field(ips_e:ipe_e,ijs_e:ije_e,ikxs:ikxe,ikys:ikye,izs:ize), ionode=0)
     ELSE
       CALL putarrnd(fidres, dset_name, field(ips_e:ipe_e,ijs_e:ije_e,ikxs:ikxe,ikys:ikye,izs:ize),  (/1,3/))
     ENDIF
     CALL attach(fidres, dset_name, 'cstep', cstep)
     CALL attach(fidres, dset_name, 'time', time)
     CALL attach(fidres, dset_name, 'jobnum', jobnum)
     CALL attach(fidres, dset_name, 'dt', dt)
     CALL attach(fidres, dset_name, 'iframe2d', iframe2d)
     CALL attach(fidres, dset_name, 'iframe5d', iframe5d)

   END SUBROUTINE write_field5d_e

   SUBROUTINE write_field5d_i(field, text)
      USE futils, ONLY: attach, putarr, putarrnd
      USE grid, ONLY: ips_i,ipe_i, ijs_i,ije_i, ikxs,ikxe, ikys,ikye, izs,ize
      USE prec_const
      IMPLICIT NONE

      COMPLEX(dp), DIMENSION(ips_i:ipe_i,ijs_i:ije_i,ikxs:ikxe,ikys:ikye,izs:ize), INTENT(IN) :: field
      CHARACTER(*), INTENT(IN) :: text

      CHARACTER(LEN=50) :: dset_name

      WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var5d", TRIM(text), iframe5d
      IF (num_procs .EQ. 1) THEN
        CALL putarr(fidres, dset_name, field(ips_i:ipe_i,ijs_i:ije_i,ikxs:ikxe,ikys:ikye,izs:ize), ionode=0)
      ELSE
        CALL putarrnd(fidres, dset_name, field(ips_i:ipe_i,ijs_i:ije_i,ikxs:ikxe,ikys:ikye,izs:ize),  (/1,3/))
      ENDIF
     CALL attach(fidres, dset_name, 'cstep', cstep)
     CALL attach(fidres, dset_name, 'time', time)
     CALL attach(fidres, dset_name, 'jobnum', jobnum)
     CALL attach(fidres, dset_name, 'dt', dt)
     CALL attach(fidres, dset_name, 'iframe2d', iframe2d)
     CALL attach(fidres, dset_name, 'iframe5d', iframe5d)

    END SUBROUTINE write_field5d_i

END SUBROUTINE diagnose_5d
