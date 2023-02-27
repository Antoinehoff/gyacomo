SUBROUTINE diagnose(kstep)
  !   Diagnostics, writing simulation state to disk
  USE basic
  USE diagnostics_par
  USE processing, ONLY: gflux_ri, hflux_xi
  IMPLICIT NONE

  INTEGER, INTENT(in) :: kstep
  CALL cpu_time(t0_diag) ! Measuring time

  !! Basic diagnose loop for reading input file, displaying advancement and ending
  IF ((kstep .EQ. 0)) THEN
    INQUIRE(unit=lu_in, name=input_fname)
    CLOSE(lu_in)
  ENDIF
  IF (kstep .GE. 0) THEN
    ! Terminal info
    IF (MOD(cstep, INT(1.0/dt)) == 0 .AND. (my_id .EQ. 0)) THEN
      ! WRITE(*,"(F6.0,A,F6.0)") time,"/",tmax
      WRITE(*,"(A,F6.0,A1,F6.0,A8,G10.2,A8,G10.2,A)")'|t/tmax = ', time,"/",tmax,'| Gxi = ',gflux_ri,'| Qxi = ',hflux_xi,'|'
    ENDIF
  ELSEIF (kstep .EQ. -1) THEN
    CALL cpu_time(finish)
     ! Display computational time cost
     IF (my_id .EQ. 0) CALL display_h_min_s(finish-start)
  END IF
  !! Specific diagnostic calls
  CALL diagnose_full(kstep)

  CALL cpu_time(t1_diag); tc_diag = tc_diag + (t1_diag - t0_diag)
END SUBROUTINE diagnose

SUBROUTINE init_outfile(comm,file0,file,fid)
  USE diagnostics_par, ONLY : write_doubleprecision, diag_par_outputinputs, input_fname
  USE basic,           ONLY : my_id, jobnum, basic_outputinputs
  USE grid,            ONLY : grid_outputinputs
  USE geometry,        ONLY : geometry_outputinputs
  USE model,           ONLY : model_outputinputs
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
  IF (my_id .EQ. 0) WRITE(*,'(3x,a,a)')  TRIM(file), ' created'
  !  basic data group
  CALL creatg(fid, "/data", "data")
  !  File group
  CALL creatg(fid, "/files", "files")
  CALL attach(fid, "/files",  "jobnum", jobnum)

  ! Add the code info and parameters to the file
  WRITE(str,'(a,i2.2)') "/data/input"
  CALL creatd(fid, 0,(/0/),TRIM(str),'Input parameters')
  CALL attach(fid, TRIM(str),     "version",  VERSION) !defined in srcinfo.h
  CALL attach(fid, TRIM(str),      "branch",   BRANCH) !defined in srcinfo.h
  CALL attach(fid, TRIM(str),      "author",   AUTHOR) !defined in srcinfo.h
  CALL attach(fid, TRIM(str),    "execdate", EXECDATE) !defined in srcinfo.h
  CALL attach(fid, TRIM(str),        "host",     HOST) !defined in srcinfo.h

  CALL basic_outputinputs(fid,str)
  CALL grid_outputinputs(fid, str)
  CALL geometry_outputinputs(fid, str)
  CALL diag_par_outputinputs(fid, str)
  CALL model_outputinputs(fid, str)
  CALL coll_outputinputs(fid, str)
  CALL initial_outputinputs(fid, str)
  CALL time_integration_outputinputs(fid, str)

  !  Save STDIN (input file) of this run
  IF(jobnum .LE. 99) THEN
     WRITE(str,'(a,i2.2)') "/files/STDIN.",jobnum
  ELSE
     WRITE(str,'(a,i3.2)') "/files/STDIN.",jobnum
  END IF
  CALL putfile(fid, TRIM(str), TRIM(input_fname),ionode=0)
END SUBROUTINE init_outfile

SUBROUTINE diagnose_full(kstep)
  USE basic
  USE grid
  USE diagnostics_par
  USE futils, ONLY: creatf, creatg, creatd, closef, putarr, putfile, attach, openf, putarrnd
  USE array
  USE model
  USE initial_par
  USE fields
  USE time_integration
  USE parallel
  USE prec_const
  USE collision, ONLY: coll_outputinputs
  USE geometry
  IMPLICIT NONE

  INTEGER, INTENT(in) :: kstep
  INTEGER, parameter  :: BUFSIZE = 2
  INTEGER :: rank = 0
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
    CALL putarr(fidres, "/data/grid/coordp_e" , parray_e_full, "p_e", ionode=0)
    CALL putarr(fidres, "/data/grid/coordj_e" , jarray_e_full, "j_e", ionode=0)
    CALL putarr(fidres, "/data/grid/coordp_i" , parray_i_full, "p_i", ionode=0)
    CALL putarr(fidres, "/data/grid/coordj_i" , jarray_i_full, "j_i", ionode=0)

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
    CALL putarrnd(fidres, "/data/metric/kernel_i",    kernel_i(ijs_i:ije_i,ikys:ikye,ikxs:ikxe,izs:ize,0:1), (/ 1, 2, 4/))

    !  var0d group (gyro transport)
    IF (nsave_0d .GT. 0) THEN
     CALL creatg(fidres, "/data/var0d", "0d profiles")
     CALL creatd(fidres, rank, dims,  "/data/var0d/time",     "Time t*c_s/R")
     CALL creatd(fidres, rank, dims, "/data/var0d/cstep", "iteration number")

     IF (write_gamma) THEN
       CALL creatd(fidres, rank, dims, "/data/var0d/gflux_ri", "Radial gyro ion transport")
       CALL creatd(fidres, rank, dims, "/data/var0d/pflux_ri", "Radial part ion transport")
       IF(KIN_E) THEN
       CALL creatd(fidres, rank, dims, "/data/var0d/gflux_re", "Radial gyro electron transport")
       CALL creatd(fidres, rank, dims, "/data/var0d/pflux_re", "Radial part electron transport")
       ENDIF
     ENDIF
     IF (write_hf) THEN
       CALL creatd(fidres, rank, dims, "/data/var0d/hflux_xi", "Radial part ion heat flux")
       IF(KIN_E) THEN
       CALL creatd(fidres, rank, dims, "/data/var0d/hflux_xe", "Radial part electron heat flux")
       ENDIF
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
      IF(KIN_E)&
      CALL creatg(fidres, "/data/var3d/Ne00", "Ne00")
      CALL creatg(fidres, "/data/var3d/Ni00", "Ni00")
      IF(KIN_E)&
      CALL creatg(fidres, "/data/var3d/Nepjz", "Nepjz")
      CALL creatg(fidres, "/data/var3d/Nipjz", "Nipjz")
     ENDIF

     IF (write_dens) THEN
       IF(KIN_E)&
      CALL creatg(fidres, "/data/var3d/dens_e", "dens_e")
      CALL creatg(fidres, "/data/var3d/dens_i", "dens_i")
     ENDIF

     IF (write_fvel) THEN
       IF(KIN_E) THEN
       CALL creatg(fidres, "/data/var3d/upar_e", "upar_e")
       CALL creatg(fidres, "/data/var3d/uper_e", "uper_e")
       ENDIF
       CALL creatg(fidres, "/data/var3d/upar_i", "upar_i")
       CALL creatg(fidres, "/data/var3d/uper_i", "uper_i")
     ENDIF

     IF (write_temp) THEN
       IF(KIN_E) THEN
       CALL creatg(fidres, "/data/var3d/Tper_e", "Tper_e")
       CALL creatg(fidres, "/data/var3d/Tpar_e", "Tpar_e")
       CALL creatg(fidres, "/data/var3d/temp_e", "temp_e")
       ENDIF
       CALL creatg(fidres, "/data/var3d/Tper_i", "Tper_i")
       CALL creatg(fidres, "/data/var3d/Tpar_i", "Tpar_i")
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
        IF(KIN_E)&
       CALL creatg(fidres, "/data/var5d/moments_e", "moments_e")
       CALL creatg(fidres, "/data/var5d/moments_i", "moments_i")
      ENDIF

      IF (write_Sapj) THEN
        IF(KIN_E)&
        CALL creatg(fidres, "/data/var5d/Sepj", "Sepj")
        CALL creatg(fidres, "/data/var5d/Sipj", "Sipj")
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
  USE model, ONLY: KIN_E

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
    CALL compute_radial_ion_transport
    CALL append(fidres, "/data/var0d/gflux_ri",gflux_ri,ionode=0)
    CALL append(fidres, "/data/var0d/pflux_ri",pflux_ri,ionode=0)
    IF(KIN_E) THEN
    CALL compute_radial_electron_transport
    CALL append(fidres, "/data/var0d/gflux_re",gflux_re,ionode=0)
    CALL append(fidres, "/data/var0d/pflux_re",pflux_re,ionode=0)
    ENDIF
  ENDIF
  IF (write_hf) THEN
    CALL compute_radial_ion_heatflux
    CALL append(fidres, "/data/var0d/hflux_xi",hflux_xi,ionode=0)
    IF(KIN_E) THEN
    CALL compute_radial_electron_heatflux
    CALL append(fidres, "/data/var0d/hflux_xe",hflux_xe,ionode=0)
    ENDIF
  ENDIF
END SUBROUTINE diagnose_0d

SUBROUTINE diagnose_3d
  USE basic
  USE futils, ONLY: append, getatt, attach, putarrnd, putarr
  USE fields
  USE array
  USE grid, ONLY: ikxs,ikxe, ikys,ikye, Nkx, Nky, local_nkx, ikx, iky, ips_e, ips_i
  USE time_integration
  USE diagnostics_par
  USE prec_const
  USE processing
  USE model, ONLY: KIN_E
  IMPLICIT NONE

  CALL append(fidres,  "/data/var3d/time",           time,ionode=0)
  CALL append(fidres, "/data/var3d/cstep", real(cstep,dp),ionode=0)
  CALL getatt(fidres,      "/data/var3d/",       "frames",iframe3d)
  iframe3d=iframe3d+1
  CALL attach(fidres,"/data/var3d/" , "frames", iframe3d)

  IF (write_phi) CALL write_field3d_kykxz(phi (ikys:ikye,ikxs:ikxe,izs:ize), 'phi')
  IF (write_phi) CALL write_field3d_kykxz(psi (ikys:ikye,ikxs:ikxe,izs:ize), 'psi')

  IF (write_Na00) THEN
    IF(KIN_E)THEN
    IF (CONTAINS_ip0_e) &
      Ne00(ikys:ikye,ikxs:ikxe,izs:ize) = moments_e(ip0_e,ij0_e,ikys:ikye,ikxs:ikxe,izs:ize,updatetlevel)
    CALL write_field3d_kykxz(Ne00(ikys:ikye,ikxs:ikxe,izs:ize), 'Ne00')
    ENDIF
    IF (CONTAINS_ip0_i) &
      Ni00(ikys:ikye,ikxs:ikxe,izs:ize) = moments_i(ip0_i,ij0_i,ikys:ikye,ikxs:ikxe,izs:ize,updatetlevel)
    CALL write_field3d_kykxz(Ni00(ikys:ikye,ikxs:ikxe,izs:ize), 'Ni00')

    CALL compute_Napjz_spectrum
    IF(KIN_E) &
    CALL write_field3d_pjz_e(Nepjz(ips_e:ipe_e,ijs_e:ije_e,izs:ize), 'Nepjz')
    CALL write_field3d_pjz_i(Nipjz(ips_i:ipe_i,ijs_i:ije_i,izs:ize), 'Nipjz')
  ENDIF

  !! Fuid moments
  IF (write_dens .OR. write_fvel .OR. write_temp) &
  CALL compute_fluid_moments

  IF (write_dens) THEN
    IF(KIN_E)&
    CALL write_field3d_kykxz(dens_e(ikys:ikye,ikxs:ikxe,izs:ize), 'dens_e')
    CALL write_field3d_kykxz(dens_i(ikys:ikye,ikxs:ikxe,izs:ize), 'dens_i')
  ENDIF

  IF (write_fvel) THEN
    IF(KIN_E)&
    CALL write_field3d_kykxz(upar_e(ikys:ikye,ikxs:ikxe,izs:ize), 'upar_e')
    CALL write_field3d_kykxz(upar_i(ikys:ikye,ikxs:ikxe,izs:ize), 'upar_i')
    IF(KIN_E)&
    CALL write_field3d_kykxz(uper_e(ikys:ikye,ikxs:ikxe,izs:ize), 'uper_e')
    CALL write_field3d_kykxz(uper_i(ikys:ikye,ikxs:ikxe,izs:ize), 'uper_i')
  ENDIF

  IF (write_temp) THEN
    IF(KIN_E)&
    CALL write_field3d_kykxz(Tpar_e(ikys:ikye,ikxs:ikxe,izs:ize), 'Tpar_e')
    CALL write_field3d_kykxz(Tpar_i(ikys:ikye,ikxs:ikxe,izs:ize), 'Tpar_i')
    IF(KIN_E)&
    CALL write_field3d_kykxz(Tper_e(ikys:ikye,ikxs:ikxe,izs:ize), 'Tper_e')
    CALL write_field3d_kykxz(Tper_i(ikys:ikye,ikxs:ikxe,izs:ize), 'Tper_i')
    IF(KIN_E)&
    CALL write_field3d_kykxz(temp_e(ikys:ikye,ikxs:ikxe,izs:ize), 'temp_e')
    CALL write_field3d_kykxz(temp_i(ikys:ikye,ikxs:ikxe,izs:ize), 'temp_i')
  ENDIF

  CONTAINS

  SUBROUTINE write_field3d_kykxz(field, text)
    USE parallel, ONLY : gather_xyz
    IMPLICIT NONE
    COMPLEX(dp), DIMENSION(ikys:ikye,ikxs:ikxe, izs:ize), INTENT(IN) :: field
    CHARACTER(*), INTENT(IN) :: text
    COMPLEX(dp), DIMENSION(1:Nky,1:Nkx,1:Nz) :: field_full
    CHARACTER(256) :: dset_name
    field_full = 0;
    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var3d", TRIM(text), iframe3d

    IF (num_procs .EQ. 1) THEN ! no data distribution
      CALL putarr(fidres, dset_name, field(ikys:ikye,ikxs:ikxe, izs:ize), ionode=0)

    ELSEIF(GATHERV_OUTPUT) THEN ! output using one node (gatherv)
      CALL gather_xyz(field(ikys:ikye,1:Nkx,izs:ize),field_full(1:Nky,1:Nkx,1:Nz))
      CALL putarr(fidres, dset_name, field_full(1:Nky,1:Nkx,1:Nz), ionode=0)
    ELSE ! output using putarrnd (very slow on marconi)
      CALL putarrnd(fidres, dset_name, field(ikys:ikye,ikxs:ikxe, izs:ize),  (/1, 1, 3/))
    ENDIF
    CALL attach(fidres, dset_name, "time", time)
  END SUBROUTINE write_field3d_kykxz

  SUBROUTINE write_field3d_pjz_i(field, text)
    USE parallel, ONLY : gather_pjz_i
    IMPLICIT NONE
    REAL(dp), DIMENSION(ips_i:ipe_i,ijs_i:ije_i,izs:ize), INTENT(IN) :: field
    REAL(dp), DIMENSION(1:Np_i,1:Nj_i,1:Nz) :: field_full
    CHARACTER(*), INTENT(IN) :: text
    CHARACTER(LEN=50) :: dset_name
    field_full = 0;
    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var3d", TRIM(text), iframe3d
    IF (num_procs .EQ. 1) THEN ! no data distribution
      CALL putarr(fidres, dset_name, field(ips_i:ipe_i,ijs_i:ije_i,izs:ize), ionode=0)
    ELSE
      CALL gather_pjz_i(field(ips_i:ipe_i,ijs_i:ije_i,izs:ize),field_full(1:Np_i,1:Nj_i,1:Nz))
      CALL putarr(fidres, dset_name, field(1:Np_i,1:Nj_i,1:Nz), ionode=0)
    ENDIF
    CALL attach(fidres, dset_name, "time", time)
  END SUBROUTINE write_field3d_pjz_i

  SUBROUTINE write_field3d_pjz_e(field, text)
    USE parallel, ONLY : gather_pjz_e
    IMPLICIT NONE
    REAL(dp), DIMENSION(ips_e:ipe_e,ijs_e:ije_e,izs:ize), INTENT(IN) :: field
    REAL(dp), DIMENSION(1:pmaxe+1,1:jmaxe+1,1:Nz) :: field_full
    CHARACTER(*), INTENT(IN) :: text
    CHARACTER(LEN=50) :: dset_name
    field_full = 0;
    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var3d", TRIM(text), iframe3d
    IF (num_procs .EQ. 1) THEN ! no data distribution
      CALL putarr(fidres, dset_name, field(ips_e:ipe_e,ijs_e:ije_e,izs:ize), ionode=0)
    ELSE
      CALL gather_pjz_e(field(ips_e:ipe_e,ijs_e:ije_e,izs:ize),field_full(1:pmaxe+1,1:jmaxe+1,1:Nz))
      CALL putarr(fidres, dset_name, field(1:Np_i,1:Nj_i,1:Nz), ionode=0)
    ENDIF
    CALL attach(fidres, dset_name, "time", time)
  END SUBROUTINE write_field3d_pjz_e

END SUBROUTINE diagnose_3d

SUBROUTINE diagnose_5d

  USE basic
  USE futils, ONLY: append, getatt, attach, putarrnd, putarr
  USE fields
  USE array!, ONLY: Sepj, Sipj
  USE grid, ONLY: ips_e,ipe_e, ips_i, ipe_i, &
                 ijs_e,ije_e, ijs_i, ije_i, &
                 Np_i, Nj_i, Np_e, Nj_e, Nky, Nkx, Nz, &
                 ikxs,ikxe,ikys,ikye,izs,ize
  USE time_integration
  USE diagnostics_par
  USE prec_const
  USE model, ONLY: KIN_E
  IMPLICIT NONE

  CALL append(fidres,  "/data/var5d/time",           time,ionode=0)
  CALL append(fidres, "/data/var5d/cstep", real(cstep,dp),ionode=0)
  CALL getatt(fidres,      "/data/var5d/",       "frames",iframe5d)
  iframe5d=iframe5d+1
  CALL attach(fidres,"/data/var5d/" , "frames", iframe5d)

  IF (write_Napj) THEN
   IF(KIN_E)&
  CALL write_field5d_e(moments_e(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,ikxs:ikxe,izs:ize,updatetlevel), 'moments_e')
  CALL write_field5d_i(moments_i(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,ikxs:ikxe,izs:ize,updatetlevel), 'moments_i')
  ENDIF

  IF (write_Sapj) THEN
   IF(KIN_E)&
   CALL write_field5d_e(Sepj(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,ikxs:ikxe,izs:ize), 'Sepj')
   CALL write_field5d_i(Sipj(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,ikxs:ikxe,izs:ize), 'Sipj')
  ENDIF

  CONTAINS

  SUBROUTINE write_field5d_e(field, text)
    USE futils, ONLY: attach, putarr, putarrnd
    USE parallel, ONLY: gather_pjxyz_e
    USE grid,   ONLY: ips_e,ipe_e, ijs_e,ije_e, ikxs,ikxe, ikys,ikye, izs,ize
    USE prec_const
    IMPLICIT NONE
    COMPLEX(dp), DIMENSION(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,ikxs:ikxe,izs:ize), INTENT(IN) :: field
    CHARACTER(*), INTENT(IN) :: text
    COMPLEX(dp), DIMENSION(1:Np_e,1:Nj_e,1:Nky,1:Nkx,1:Nz) :: field_full
    CHARACTER(LEN=50) :: dset_name
    field_full = 0;
    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var5d", TRIM(text), iframe5d
    IF (num_procs .EQ. 1) THEN
     CALL putarr(fidres, dset_name, field(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,ikxs:ikxe,izs:ize), ionode=0)
   ELSEIF(GATHERV_OUTPUT) THEN ! output using one node (gatherv)
     CALL gather_pjxyz_e(field(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,ikxs:ikxe,izs:ize),&
                         field_full(1:Np_e,1:Nj_e,1:Nky,1:Nkx,1:Nz))
     CALL putarr(fidres, dset_name, field_full(1:Np_i,1:Nj_i,1:Nky,1:Nkx,1:Nz), ionode=0)
   ELSE
     CALL putarrnd(fidres, dset_name, field(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,ikxs:ikxe,izs:ize),  (/1,3,5/))
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
    USE parallel, ONLY: gather_pjxyz_i
    USE grid, ONLY: ips_i,ipe_i, ijs_i,ije_i, ikxs,ikxe, ikys,ikye, izs,ize
    USE prec_const
    IMPLICIT NONE
    COMPLEX(dp), DIMENSION(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,ikxs:ikxe,izs:ize), INTENT(IN) :: field
    CHARACTER(*), INTENT(IN) :: text
    COMPLEX(dp), DIMENSION(1:Np_i,1:Nj_i,1:Nky,1:Nkx,1:Nz) :: field_full
    CHARACTER(LEN=50) :: dset_name
    field_full = 0;
    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var5d", TRIM(text), iframe5d
    IF (num_procs .EQ. 1) THEN
      CALL putarr(fidres, dset_name, field(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,ikxs:ikxe,izs:ize), ionode=0)
    ELSEIF(GATHERV_OUTPUT) THEN ! output using one node (gatherv)
      CALL gather_pjxyz_i(field(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,ikxs:ikxe,izs:ize),&
                        field_full(1:Np_i,1:Nj_i,1:Nky,1:Nkx,1:Nz))
      CALL putarr(fidres, dset_name, field_full(1:Np_i,1:Nj_i,1:Nky,1:Nkx,1:Nz), ionode=0)
    ELSE
      CALL putarrnd(fidres, dset_name, field(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,ikxs:ikxe,izs:ize),  (/1,3,5/))
    ENDIF
    CALL attach(fidres, dset_name, 'cstep', cstep)
    CALL attach(fidres, dset_name, 'time', time)
    CALL attach(fidres, dset_name, 'jobnum', jobnum)
    CALL attach(fidres, dset_name, 'dt', dt)
    CALL attach(fidres, dset_name, 'iframe2d', iframe2d)
    CALL attach(fidres, dset_name, 'iframe5d', iframe5d)

  END SUBROUTINE write_field5d_i
END SUBROUTINE diagnose_5d

SUBROUTINE spit_snapshot_check
  USE fields, ONLY: phi
  USE grid, ONLY: ikxs,ikxe,Nkx,ikys,ikye,Nky,izs,ize,Nz
  USE parallel, ONLY: gather_xyz
  USE basic
  IMPLICIT NONE
  LOGICAL :: file_exist
  INTEGER :: fid_check, ikx, iky, iz
  CHARACTER(256) :: check_filename
  COMPLEX(dp), DIMENSION(1:Nky,1:Nkx,1:Nz) :: field_to_check
  !! Spit a snapshot of PHI if requested (triggered by creating a file named "check_phi")
  INQUIRE(file='check_phi', exist=file_exist)
  IF( file_exist ) THEN
     IF(my_id.EQ. 0) WRITE(*,*) 'Check file found -> gather phi..'
     CALL gather_xyz(phi(ikys:ikye,ikxs:ikxe,izs:ize), field_to_check)
     IF(my_id.EQ. 0) THEN
       WRITE(check_filename,'(a16)') 'check_phi.out'
       OPEN(fid_check, file=check_filename, form='formatted')
       WRITE(*,*) 'Check file found -> output phi ..'
       WRITE(fid_check,*) Nky, Nkx, Nz
       DO iky = 1,Nky; DO ikx = 1, Nkx; DO iz = 1,Nz
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
