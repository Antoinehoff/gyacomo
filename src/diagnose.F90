SUBROUTINE diagnose(kstep)
  !   Diagnostics, writing simulation state to disk
  USE basic, ONLY: t0_diag,t1_diag,tc_diag, lu_in, finish, start, cstep, dt, time, tmax, display_h_min_s
  USE diagnostics_par, ONLY: input_fname
  USE processing, ONLY: pflux_x, hflux_x
  USE parallel,   ONLY: my_id
  IMPLICIT NONE

  INTEGER, INTENT(in) :: kstep
  CALL cpu_time(t0_diag) ! Measuring time

  !! Basic diagnose loop for reading input file, displaying advancement and ending
  IF ((kstep .EQ. 0)) THEN
    INQUIRE(unit=lu_in, name=input_fname)
    CLOSE(lu_in)
  ENDIF
  !! End diag
  IF (kstep .EQ. -1) THEN
    CALL cpu_time(finish)
     ! Display computational time cost
     CALL display_h_min_s(finish-start)
     ! Show last state transport values
     IF (my_id .EQ. 0) &
      WRITE(*,"(A,G10.2,A8,G10.2,A)") 'Final transport values : | Gxi = ',pflux_x(1),'| Qxi = ',hflux_x(1),'|'
  END IF
  !! Specific diagnostic calls
  CALL diagnose_full(kstep)
  ! Terminal info
  IF ((kstep .GT. 0) .AND. (MOD(cstep, INT(1.0/dt)) == 0) .AND. (my_id .EQ. 0)) THEN
    WRITE(*,"(A,F6.0,A1,F6.0,A8,G10.2,A8,G10.2,A)")'|t/tmax = ', time,"/",tmax,'| Gxi = ',pflux_x(1),'| Qxi = ',hflux_x(1),'|'
  ENDIF
  CALL cpu_time(t1_diag); tc_diag = tc_diag + (t1_diag - t0_diag)
END SUBROUTINE diagnose

SUBROUTINE init_outfile(comm,file0,file,fid)
  USE diagnostics_par, ONLY : write_doubleprecision, diag_par_outputinputs, input_fname
  USE basic,           ONLY : speak, jobnum, basic_outputinputs
  USE grid,            ONLY : grid_outputinputs
  USE geometry,        ONLY : geometry_outputinputs
  USE model,           ONLY : model_outputinputs
  USE species,         ONLY : species_outputinputs
  USE collision,       ONLY : coll_outputinputs
  USE initial_par,     ONLY : initial_outputinputs
  USE time_integration,ONLY : time_integration_outputinputs
  USE futils,          ONLY : creatf, creatg, creatd, attach, putfile
  IMPLICIT NONE
  !input
  INTEGER,            INTENT(IN)    :: comm
  CHARACTER(len=256), INTENT(IN)    :: file0
  CHARACTER(len=256), INTENT(OUT)   :: file
  INTEGER,            INTENT(OUT)   :: fid
  CHARACTER(len=256)                :: str
  INCLUDE 'srcinfo.h'

  ! Writing output filename
  WRITE(file,'(a,a1,i2.2,a3)') TRIM(file0)   ,'_',jobnum,'.h5'
  !                      1.1   Initial run
  ! Main output file creation
  IF (write_doubleprecision) THEN
    CALL creatf(file, fid, real_prec='d', mpicomm=comm)
  ELSE
    CALL creatf(file, fid, mpicomm=comm)
  END IF
  CALL speak(TRIM(file)//' created')
  !  basic data group
  CALL creatg(fid, "/data", "data")
  !  File group
  CALL creatg(fid, "/files", "files")
  CALL attach(fid, "/files",  "jobnum", jobnum)

  ! Add the code info and parameters to the file
  CALL creatg(fid, "/data/input", "input")
  CALL creatd(fid, 0,(/0/),"/data/input/codeinfo",'Code Information')
  CALL attach(fid, TRIM(str),     "version",  VERSION) !defined in srcinfo.h
  CALL attach(fid, TRIM(str),      "branch",   BRANCH) !defined in srcinfo.h
  CALL attach(fid, TRIM(str),      "author",   AUTHOR) !defined in srcinfo.h
  CALL attach(fid, TRIM(str),    "execdate", EXECDATE) !defined in srcinfo.h
  CALL attach(fid, TRIM(str),        "host",     HOST) !defined in srcinfo.h
  CALL basic_outputinputs(fid)
  CALL grid_outputinputs(fid)
  CALL geometry_outputinputs(fid)
  CALL diag_par_outputinputs(fid)
  CALL model_outputinputs(fid)
  CALL species_outputinputs(fid)
  CALL coll_outputinputs(fid)
  CALL initial_outputinputs(fid)
  CALL time_integration_outputinputs(fid)

  !  Save STDIN (input file) of this run
  IF(jobnum .LE. 99) THEN
     WRITE(str,'(a,i2.2)') "/files/STDIN.",jobnum
  ELSE
     WRITE(str,'(a,i3.2)') "/files/STDIN.",jobnum
  END IF
  CALL putfile(fid, TRIM(str), TRIM(input_fname),ionode=0)
END SUBROUTINE init_outfile

SUBROUTINE diagnose_full(kstep)
  USE basic,           ONLY: speak,&
                             cstep,iframe0d,iframe2d,iframe3d,iframe5d,&
                             start,finish,crashed
  USE grid,            ONLY: &
    ias,iae, &
    parray_full,pmax,&
    ijs,ije,jarray_full,jmax,&
    ikys,ikye,kyarray_full,&
    ikxs,ikxe,kxarray_full,&
    izs,ize,zarray_full
  USE diagnostics_par
  USE futils,          ONLY: creatf, creatg, creatd, closef, putarr, putfile, attach, openf, putarrnd
  USE array
  USE model,           ONLY:
  USE initial_par,     ONLY:
  USE fields,          ONLY:
  USE time_integration,ONLY:
  USE parallel,        ONLY: my_id, comm0
  USE prec_const,      ONLY:
  USE collision,       ONLY: coll_outputinputs
  USE geometry,        ONLY: gxx,gxy,gyy,gxz,gyz,gzz,hatR,hatZ,hatB,dBdx,dBdy,dBdz,Jacobian,gradz_coeff,Ckxky
  IMPLICIT NONE

  INTEGER, INTENT(in) :: kstep
  INTEGER, parameter  :: BUFSIZE = 2
  INTEGER :: rank = 0, ierr
  INTEGER :: dims(1) = (/0/)
  !____________________________________________________________________________
  !                   1.   Initial diagnostics

  IF ((kstep .EQ. 0)) THEN
    CALL init_outfile(comm0,   resfile0,resfile,fidres)

    ! Profiler time measurement
    CALL creatg(fidres, "/profiler", "performance analysis")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_rhs",        "cumulative rhs computation time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_poisson",    "cumulative poisson computation time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_Sapj",       "cumulative Sapj computation time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_coll",       "cumulative collision computation time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_process",    "cumulative process computation time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_adv_field",  "cumulative adv. fields computation time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_ghost",       "cumulative communication time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_clos",       "cumulative closure computation time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_checkfield", "cumulative checkfield computation time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_diag",       "cumulative sym computation time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_step",       "cumulative total step computation time")
    CALL creatd(fidres, 0, dims, "/profiler/time",          "current simulation time")

    ! Grid info
    CALL creatg(fidres, "/data/grid", "Grid data")
    CALL putarr(fidres, "/data/grid/coordkx",   kxarray_full,  "kx*rho_s0", ionode=0)
    CALL putarr(fidres, "/data/grid/coordky",   kyarray_full,  "ky*rho_s0", ionode=0)
    CALL putarr(fidres, "/data/grid/coordz",    zarray_full,   "z/R", ionode=0)
    CALL putarr(fidres, "/data/grid/coordp" ,   parray_full,   "p", ionode=0)
    CALL putarr(fidres, "/data/grid/coordj" ,   jarray_full,   "j", ionode=0)

    ! Metric info
    CALL   creatg(fidres, "/data/metric", "Metric data")
    CALL putarrnd(fidres, "/data/metric/gxx",            gxx(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidres, "/data/metric/gxy",            gxy(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidres, "/data/metric/gxz",            gxz(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidres, "/data/metric/gyy",            gyy(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidres, "/data/metric/gyz",            gyz(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidres, "/data/metric/gzz",            gzz(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidres, "/data/metric/hatR",          hatR(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidres, "/data/metric/hatZ",          hatZ(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidres, "/data/metric/hatB",          hatB(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidres, "/data/metric/dBdx",      dBdx(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidres, "/data/metric/dBdy",      dBdy(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidres, "/data/metric/dBdz",      dBdz(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidres, "/data/metric/Jacobian",    Jacobian(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidres, "/data/metric/gradz_coeff", gradz_coeff(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidres, "/data/metric/Ckxky",       Ckxky(ikys:ikye,ikxs:ikxe,izs:ize,0:1), (/1, 1, 3/))
    CALL putarrnd(fidres, "/data/metric/kernel",    kernel(ias:iae,ijs:ije,ikys:ikye,ikxs:ikxe,izs:ize,0:1), (/1, 1, 2, 4/))

    !  var0d group (gyro transport)
    IF (nsave_0d .GT. 0) THEN
     CALL creatg(fidres, "/data/var0d", "0d profiles")
     CALL creatd(fidres, rank, dims,  "/data/var0d/time",     "Time t*c_s/R")
     CALL creatd(fidres, rank, dims, "/data/var0d/cstep", "iteration number")

     IF (write_gamma) THEN
       CALL creatd(fidres, rank, dims, "/data/var0d/gflux_x", "Radial gyro transport")
       CALL creatd(fidres, rank, dims, "/data/var0d/pflux_x", "Radial part transport")
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
     IF (write_phi) CALL creatg(fidres, "/data/var3d/psi", "psi")

     IF (write_Na00) THEN
      CALL creatg(fidres, "/data/var3d/Na00", "Na00")
      CALL creatg(fidres, "/data/var3d/Nepjz", "Napjz")
     ENDIF

     IF (write_dens) THEN
      CALL creatg(fidres, "/data/var3d/dens", "dens_e")
     ENDIF

     IF (write_fvel) THEN
       CALL creatg(fidres, "/data/var3d/upar", "upar")
       CALL creatg(fidres, "/data/var3d/uper", "uper")
     ENDIF

     IF (write_temp) THEN
       CALL creatg(fidres, "/data/var3d/Tper_e", "Tper")
       CALL creatg(fidres, "/data/var3d/Tpar_e", "Tpar")
       CALL creatg(fidres, "/data/var3d/temp_e", "temp")
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
       CALL creatg(fidres, "/data/var5d/moments", "moments")
      ENDIF

      IF (write_Sapj) THEN
        CALL creatg(fidres, "/data/var5d/Sapj", "Sapj")
      ENDIF

      IF (cstep==0) THEN
       iframe5d=0
      END IF
      CALL attach(fidres,"/data/var5d/" , "frames", iframe5d)
    END IF
  ENDIF

  !_____________________________________________________________________________
  !                   2.   Periodic diagnostics
  !
  IF (kstep .GE. 0) THEN

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


     !                       2.3   3d profiles
     IF (nsave_3d .GT. 0) THEN
        IF (MOD(cstep, nsave_3d) == 0) THEN
          CALL diagnose_3d
          ! Looks at the folder if the file check_phi exists and spits a snapshot
          ! of the current electrostatic potential in a basic text file
          CALL spit_snapshot_check
        ENDIF
     ENDIF

     !                       2.4   5d profiles
     IF (nsave_5d .GT. 0 .AND. cstep .GT. 0) THEN
        IF (MOD(cstep, nsave_5d) == 0) THEN
           CALL diagnose_5d
        END IF
     END IF

  !_____________________________________________________________________________
  !                   3.   Final diagnostics

  ELSEIF (kstep .EQ. -1) THEN
     CALL attach(fidres, "/data/input","cpu_time",finish-start)

     ! make a checkpoint at last timestep if not crashed
     IF(.NOT. crashed) THEN
       IF(my_id .EQ. 0) write(*,*) 'Saving last state'
       IF (nsave_5d .GT. 0) &
       CALL diagnose_5d
     ENDIF

     !   Close all diagnostic files
     CALL mpi_barrier(MPI_COMM_WORLD, ierr)
     CALL closef(fidres)

  END IF
END SUBROUTINE diagnose_full

!!-------------- Auxiliary routines -----------------!!
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
  CALL append(fidres, "/profiler/Tc_process",      tc_process,ionode=0)
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
    CALL compute_radial_transport
    CALL append(fidres, "/data/var0d/gflux_x",gflux_x(1),ionode=0)
    CALL append(fidres, "/data/var0d/pflux_x",pflux_x(1),ionode=0)
  ENDIF
  IF (write_hf) THEN
    CALL compute_radial_heatflux
    CALL append(fidres, "/data/var0d/hflux_x",hflux_x(1),ionode=0)
  ENDIF
END SUBROUTINE diagnose_0d

SUBROUTINE diagnose_3d
  USE basic
  USE futils, ONLY: append, getatt, attach, putarrnd, putarr
  USE fields, ONLY: phi, psi, moments
  USE array,  ONLY: Na00,Napjz,dens,upar,uper,Tpar,Tper,temp
  USE grid, ONLY: CONTAINSp0, ip0,ij0, &
                  total_np, total_nj, total_nky, total_nkx, total_nz, &
                  local_np, local_nj, local_nky, local_nkx, local_nz, &
                  ngz
  USE time_integration, ONLY: updatetlevel
  USE diagnostics_par
  USE prec_const
  USE processing, ONLY: compute_fluid_moments, compute_Napjz_spectrum
  IMPLICIT NONE

  CALL append(fidres,  "/data/var3d/time",           time,ionode=0)
  CALL append(fidres, "/data/var3d/cstep", real(cstep,dp),ionode=0)
  CALL getatt(fidres,      "/data/var3d/",       "frames",iframe3d)
  iframe3d=iframe3d+1
  CALL attach(fidres,"/data/var3d/" , "frames", iframe3d)

  IF (write_phi) CALL write_field3d_kykxz(phi (:,:,1+ngz/2:local_nz+ngz/2), 'phi')
  IF (write_phi) CALL write_field3d_kykxz(psi (:,:,1+ngz/2:local_nz+ngz/2), 'psi')

  IF (write_Na00) THEN
    IF (CONTAINSp0) &
      Na00(1,:,:,:) = moments(1,ip0,ij0,:,:,1+ngz/2:local_nz+ngz/2,updatetlevel)
    CALL write_field3d_kykxz(Na00(1,:,:,:), 'Na00')

    CALL compute_Napjz_spectrum
    CALL write_field3d_pjz(Napjz(1,:,:,:), 'Napjz')
  ENDIF

  !! Fuid moments
  IF (write_dens .OR. write_fvel .OR. write_temp) &
  CALL compute_fluid_moments

  IF (write_dens) THEN
    CALL write_field3d_kykxz(dens(1,:,:,:), 'dens')
  ENDIF

  IF (write_fvel) THEN
    CALL write_field3d_kykxz(upar(1,:,:,:), 'upar')
    CALL write_field3d_kykxz(uper(1,:,:,:), 'uper')
  ENDIF

  IF (write_temp) THEN
    CALL write_field3d_kykxz(Tpar(1,:,:,:), 'Tpar')
    CALL write_field3d_kykxz(Tper(1,:,:,:), 'Tper')
    CALL write_field3d_kykxz(temp(1,:,:,:), 'temp')
  ENDIF

  CONTAINS
    SUBROUTINE write_field3d_kykxz(field, text)
      USE basic,    ONLY : GATHERV_OUTPUT
      USE parallel, ONLY : gather_xyz, num_procs
      IMPLICIT NONE
      COMPLEX(dp), DIMENSION(local_nky,local_nkx,local_nz), INTENT(IN) :: field
      CHARACTER(*), INTENT(IN) :: text
      COMPLEX(dp), DIMENSION(total_nky,total_nkx,total_nz) :: field_full
      CHARACTER(256) :: dset_name
      field_full = 0;
      WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var3d", TRIM(text), iframe3d

      IF (num_procs .EQ. 1) THEN ! no data distribution
        CALL putarr(fidres, dset_name, field, ionode=0)

      ELSEIF(GATHERV_OUTPUT) THEN ! output using one node (gatherv)
        CALL gather_xyz(field,field_full,local_nky,total_nky,total_nkx,local_nz,total_nz)
        CALL putarr(fidres, dset_name, field_full, ionode=0)
      ELSE ! output using putarrnd (very slow on marconi)
        CALL putarrnd(fidres, dset_name, field,  (/1, 1, 3/))
      ENDIF
      CALL attach(fidres, dset_name, "time", time)
    END SUBROUTINE write_field3d_kykxz

    SUBROUTINE write_field3d_pjz(field, text)
      USE parallel, ONLY : gather_pjz, num_procs
      IMPLICIT NONE
      REAL(dp), DIMENSION(local_np,local_nj,local_nz), INTENT(IN) :: field
      REAL(dp), DIMENSION(total_np,total_nj,total_nz) :: field_full
      CHARACTER(*), INTENT(IN) :: text
      CHARACTER(LEN=50) :: dset_name
      field_full = 0;
      WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var3d", TRIM(text), iframe3d
      IF (num_procs .EQ. 1) THEN ! no data distribution
        CALL putarr(fidres, dset_name, field, ionode=0)
      ELSE
        CALL gather_pjz(field,field_full,local_np,total_np,total_nj,local_nz,total_nz)
        CALL putarr(fidres, dset_name, field_full, ionode=0)
      ENDIF
      CALL attach(fidres, dset_name, "time", time)
    END SUBROUTINE write_field3d_pjz

END SUBROUTINE diagnose_3d

SUBROUTINE diagnose_5d

  USE basic,  ONLY: time, iframe5d,cstep
  USE futils, ONLY: append, getatt, attach, putarrnd, putarr
  USE fields, ONLY: moments
  USE array,  ONLY: Sapj
  USE grid,   ONLY:ips, ipe, ijs, ije, &
                   total_np, total_nj, total_nky, total_nkx, total_nz, &
                   local_np, local_nj, local_nky, local_nkx, local_nz, &
                   ngp, ngj, ngz,&
                 ikxs,ikxe,ikys,ikye,izs,ize
  USE time_integration, ONLY: updatetlevel
  USE diagnostics_par
  USE prec_const, ONLY: dp
  IMPLICIT NONE

  CALL append(fidres,  "/data/var5d/time",           time,ionode=0)
  CALL append(fidres, "/data/var5d/cstep", real(cstep,dp),ionode=0)
  CALL getatt(fidres,      "/data/var5d/",       "frames",iframe5d)
  iframe5d=iframe5d+1
  CALL attach(fidres,"/data/var5d/" , "frames", iframe5d)

  IF (write_Napj) THEN
  CALL write_field5d(moments(1,1+ngp/2:local_np+ngp/2,1+ngj/2:local_nj+ngj/2,&
                           :,:,1+ngz/2:local_nz+ngz/2,updatetlevel), 'moments')
  ENDIF

  IF (write_Sapj) THEN
   CALL write_field5d(Sapj(1,:,:,:,:,:), 'Sapj')
  ENDIF

  CONTAINS

  SUBROUTINE write_field5d(field, text)
    USE basic,    ONLY: GATHERV_OUTPUT, jobnum, dt
    USE futils,   ONLY: attach, putarr, putarrnd
    USE parallel, ONLY: gather_pjxyz, num_procs
    USE prec_const, ONLY: dp
    IMPLICIT NONE
    COMPLEX(dp), DIMENSION(local_np,local_nj,local_nky,local_nkx,local_nz), INTENT(IN) :: field
    CHARACTER(*), INTENT(IN) :: text
    COMPLEX(dp), DIMENSION(total_np,total_nj,total_nky,total_nkx,total_nz) :: field_full
    CHARACTER(LEN=50) :: dset_name
    field_full = 0;
    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var5d", TRIM(text), iframe5d
    IF (num_procs .EQ. 1) THEN
      CALL putarr(fidres, dset_name, field(ips:ipe,ijs:ije,ikys:ikye,ikxs:ikxe,izs:ize), ionode=0)
    ELSEIF(GATHERV_OUTPUT) THEN ! output using one node (gatherv)
      CALL gather_pjxyz(field,field_full,local_np,total_np,total_nj,local_nky,total_nky,total_nkx,local_nz,total_nz)
      CALL putarr(fidres, dset_name, field_full(1:total_np,1:total_nj,1:total_nky,1:total_nkx,1:total_nz), ionode=0)
    ELSE
      CALL putarrnd(fidres, dset_name, field(ips:ipe,ijs:ije,ikys:ikye,ikxs:ikxe,izs:ize),  (/1,3,5/))
    ENDIF
    CALL attach(fidres, dset_name, 'cstep', cstep)
    CALL attach(fidres, dset_name, 'time', time)
    CALL attach(fidres, dset_name, 'jobnum', jobnum)
    CALL attach(fidres, dset_name, 'dt', dt)
    CALL attach(fidres, dset_name, 'iframe5d', iframe5d)
  END SUBROUTINE write_field5d
END SUBROUTINE diagnose_5d

SUBROUTINE spit_snapshot_check
  USE fields, ONLY: phi
  USE grid, ONLY: total_nkx,total_nky,total_nz,&
                  local_nky,local_nz, ngz
  USE parallel, ONLY: gather_xyz, my_id
  USE basic
  USE prec_const, ONLY: dp
  IMPLICIT NONE
  LOGICAL :: file_exist
  INTEGER :: fid_check, ikx, iky, iz
  CHARACTER(256) :: check_filename
  COMPLEX(dp), DIMENSION(total_nky,total_nkx,total_nz) :: field_to_check
  !! Spit a snapshot of PHI if requested (triggered by creating a file named "check_phi")
  INQUIRE(file='check_phi', exist=file_exist)
  IF( file_exist ) THEN
     IF(my_id.EQ. 0) WRITE(*,*) 'Check file found -> gather phi..'
     CALL gather_xyz(phi(:,:,1+Ngz/2:local_nz+Ngz/2), field_to_check,local_nky,total_nky,total_nkx,local_nz,total_nz)
     IF(my_id.EQ. 0) THEN
       WRITE(check_filename,'(a16)') 'check_phi.out'
       OPEN(fid_check, file=check_filename, form='formatted')
       WRITE(*,*) 'Check file found -> output phi ..'
       WRITE(fid_check,*) total_nky, total_nkx, total_nz
       DO iky = 1,total_nky; DO ikx = 1, total_nkx; DO iz = 1,total_nz
         WRITE(fid_check,*) real(field_to_check(iky,ikx,iz)), ',' , imag(field_to_check(iky,ikx,iz))
       ENDDO; ENDDO; ENDDO
       CLOSE(fid_check)
       WRITE(*,*) 'Check file found -> done.'
       ! delete the check_phi flagfile
       OPEN(fid_check, file='check_phi')
       CLOSE(fid_check, status='delete')
     ENDIF
  ENDIF
END SUBROUTINE spit_snapshot_check
