SUBROUTINE diagnose(kstep)
  !   Diagnostics, writing simulation state to disk
  USE basic,           ONLY: lu_in, chrono_runt, cstep, dt, time, tmax, display_h_min_s
  USE diagnostics_par, ONLY: input_fname, diag_mode
  USE processing,      ONLY: pflux_x, hflux_x
  USE parallel,        ONLY: my_id
  IMPLICIT NONE
  INTEGER, INTENT(in) :: kstep
  !! Basic diagnose loop for reading input file, displaying advancement and ending
  IF ((kstep .EQ. 0)) THEN
    INQUIRE(unit=lu_in, name=input_fname)
    CLOSE(lu_in)
  ENDIF
  !! End diag
  IF (kstep .EQ. -1) THEN
     ! Display total run time
     CALL display_h_min_s(chrono_runt%ttot)
  END IF
  !! Specific diagnostic calls
  SELECT CASE(diag_mode)
  CASE('full')
    CALL diagnose_full(kstep)
  CASE('txtonly')
    CALL diagnose_txtonly(kstep)
  END SELECT

  ! Terminal info
  IF ((kstep .GE. 0) .AND. (MOD(cstep, INT(1.0/dt)) == 0) .AND. (my_id .EQ. 0)) THEN
    WRITE(*,"(A,F6.0,A1,F6.0,A8,G10.2,A8,G10.2,A)")'|t/tmax = ', time,"/",tmax,'| Gxi = ',pflux_x(1),'| Qxi = ',hflux_x(1),'|'
  ENDIF
END SUBROUTINE diagnose

SUBROUTINE init_outfile(comm,file0,file,fid)
  USE diagnostics_par, ONLY: write_doubleprecision, diag_par_outputinputs, input_fname
  USE basic,           ONLY: speak, jobnum, basic_outputinputs
  USE grid,            ONLY: grid_outputinputs
  USE geometry,        ONLY: geometry_outputinputs
  USE model,           ONLY: model_outputinputs
  USE closure,         ONLY: closure_outputinputs
  USE species,         ONLY: species_outputinputs
  USE collision,       ONLY: coll_outputinputs
  USE initial_par,     ONLY: initial_outputinputs
  USE time_integration,ONLY: time_integration_outputinputs
  USE futils,          ONLY: creatf, creatg, creatd, attach, putfile
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
  CALL attach(fid, "/data/input/codeinfo",  "version",  VERSION) !defined in srcinfo.h
  CALL attach(fid, "/data/input/codeinfo",   "branch",   BRANCH) !defined in srcinfo.h
  CALL attach(fid, "/data/input/codeinfo",   "author",   AUTHOR) !defined in srcinfo.h
  CALL attach(fid, "/data/input/codeinfo", "execdate", EXECDATE) !defined in srcinfo.h
  CALL attach(fid, "/data/input/codeinfo",     "host",     HOST) !defined in srcinfo.h
  CALL            basic_outputinputs(fid)
  CALL             grid_outputinputs(fid)
  CALL         geometry_outputinputs(fid)
  CALL         diag_par_outputinputs(fid)
  CALL            model_outputinputs(fid)
  CALL          closure_outputinputs(fid)
  CALL          species_outputinputs(fid)
  CALL             coll_outputinputs(fid)
  CALL          initial_outputinputs(fid)
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
  USE basic,           ONLY: speak,chrono_runt,&
                             cstep,iframe0d,iframe3d,iframe5d,crashed
  USE grid,            ONLY: &
    parray_full,pmax,jarray_full,jmax, kparray, &
    kyarray_full,kxarray_full,zarray_full, ngz, total_nz, local_nz, ieven,&
    local_Nky, total_nky, local_nkx, total_nkx
  USE geometry, ONLY: gxx, gxy, gxz, gyy, gyz, gzz, &
                      hatR, hatZ, hatB, dBdx, dBdy, dBdz, Jacobian, gradz_coeff
  USE diagnostics_par
  USE futils,          ONLY: creatf, creatg, creatd, closef, putarr, putfile, attach, openf!, putarrnd ! Routine de merde, jamais l'utiliser
  USE array
  USE model,           ONLY: EM
  USE parallel,        ONLY: my_id, comm0, gather_z, gather_xyz_real
  USE collision,       ONLY: coll_outputinputs
  IMPLICIT NONE
  INTEGER, INTENT(in) :: kstep
  INTEGER, parameter  :: BUFSIZE = 2
  REAL(xp), DIMENSION(total_nz) :: Az_full ! full z array for metric output
  REAL(xp), DIMENSION(total_nky,total_nkx,total_nz) :: kp_full
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
    CALL creatd(fidres, 0, dims, "/profiler/Tc_grad",       "cumulative grad computation time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_nadiab",     "cumulative nadiab moments computation time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_adv_field",  "cumulative adv. fields computation time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_ghost",       "cumulative communication time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_clos",       "cumulative closure computation time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_checkfield", "cumulative checkfield computation time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_diag",       "cumulative sym computation time")
    CALL creatd(fidres, 0, dims, "/profiler/Tc_step",       "cumulative total step computation time")
    CALL creatd(fidres, 0, dims, "/profiler/time",          "current simulation time")
#ifdef TEST_SVD
    CALL creatd(fidres, 0, (/0/), "/profiler/Tc_CLA", "cumulative total CLA computation time")
#endif
    ! Grid info
    CALL creatg(fidres, "/data/grid", "Grid data")
    CALL putarr(fidres, "/data/grid/coordkx",   kxarray_full,  "kx*rho_s0", ionode=0)
    CALL putarr(fidres, "/data/grid/coordky",   kyarray_full,  "ky*rho_s0", ionode=0)
    CALL putarr(fidres, "/data/grid/coordz",    zarray_full,   "z/R", ionode=0)
    CALL putarr(fidres, "/data/grid/coordp" ,   parray_full,   "p", ionode=0)
    CALL putarr(fidres, "/data/grid/coordj" ,   jarray_full,   "j", ionode=0)
    CALL gather_xyz_real(kparray(1:local_Nky,1:local_Nkx,1:local_nz,ieven),kp_full,local_nky,total_nky,total_nkx,local_nz,total_nz)
    CALL putarr(fidres, "/data/grid/coordkp" ,      kp_full,   "kp", ionode=0)
    ! Metric info
    CALL   creatg(fidres, "/data/metric", "Metric data")
    CALL gather_z(gxx((1+ngz/2):(local_nz+ngz/2),ieven),Az_full,local_nz,total_nz)
    CALL putarr(fidres, "/data/metric/gxx", Az_full, "gxx", ionode =0)
    CALL gather_z(gxy((1+ngz/2):(local_nz+ngz/2),ieven),Az_full,local_nz,total_nz)
    CALL putarr(fidres, "/data/metric/gxy", Az_full, "gxy", ionode =0)
    CALL gather_z(gxz((1+ngz/2):(local_nz+ngz/2),ieven),Az_full,local_nz,total_nz)
    CALL putarr(fidres, "/data/metric/gxz", Az_full, "gxz", ionode =0)
    CALL gather_z(gyy((1+ngz/2):(local_nz+ngz/2),ieven),Az_full,local_nz,total_nz)
    CALL putarr(fidres, "/data/metric/gyy", Az_full, "gyy", ionode =0)
    CALL gather_z(gyz((1+ngz/2):(local_nz+ngz/2),ieven),Az_full,local_nz,total_nz)
    CALL putarr(fidres, "/data/metric/gyz", Az_full, "gyz", ionode =0)
    CALL gather_z(gzz((1+ngz/2):(local_nz+ngz/2),ieven),Az_full,local_nz,total_nz)
    CALL putarr(fidres, "/data/metric/gzz", Az_full, "gzz", ionode =0)
    CALL gather_z(hatR((1+ngz/2):(local_nz+ngz/2),ieven),Az_full,local_nz,total_nz)
    CALL putarr(fidres, "/data/metric/hatR", Az_full, "hatR", ionode =0)
    CALL gather_z(hatZ((1+ngz/2):(local_nz+ngz/2),ieven),Az_full,local_nz,total_nz)
    CALL putarr(fidres, "/data/metric/hatZ", Az_full, "hatZ", ionode =0)
    CALL gather_z(hatB((1+ngz/2):(local_nz+ngz/2),ieven),Az_full,local_nz,total_nz)
    CALL putarr(fidres, "/data/metric/hatB", Az_full, "hatB", ionode =0)
    CALL gather_z(dBdx((1+ngz/2):(local_nz+ngz/2),ieven),Az_full,local_nz,total_nz)
    CALL putarr(fidres, "/data/metric/dBdx", Az_full, "dBdx", ionode =0)
    CALL gather_z(dBdy((1+ngz/2):(local_nz+ngz/2),ieven),Az_full,local_nz,total_nz)
    CALL putarr(fidres, "/data/metric/dBdy", Az_full, "dBdy", ionode =0)
    CALL gather_z(dBdz((1+ngz/2):(local_nz+ngz/2),ieven),Az_full,local_nz,total_nz)
    CALL putarr(fidres, "/data/metric/dBdz", Az_full, "dBdz", ionode =0)
    CALL gather_z(Jacobian((1+ngz/2):(local_nz+ngz/2),ieven),Az_full,local_nz,total_nz)
    CALL putarr(fidres, "/data/metric/Jacobian", Az_full, "Jacobian", ionode =0)
    CALL gather_z(gradz_coeff((1+ngz/2):(local_nz+ngz/2),ieven),Az_full,local_nz,total_nz)
    CALL putarr(fidres, "/data/metric/gradz_coeff", Az_full, "gradz_coeff", ionode =0)
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
      CALL creatd(fidres, rank, dims, "/data/var0d/hflux_x", "Radial part heat flux")
     ENDIF
     IF (cstep==0) THEN
       iframe0d=0
     ENDIF
     CALL attach(fidres,"/data/var0d/" , "frames", iframe0d)
    END IF
    !  var2d group (gyro transport)
    IF (nsave_0d .GT. 0) THEN
      CALL creatg(fidres, "/data/var2d", "2d profiles")
      CALL creatd(fidres, rank, dims,  "/data/var2d/time",     "Time t*c_s/R")
      CALL creatd(fidres, rank, dims, "/data/var2d/cstep", "iteration number")
#ifdef TEST_SVD
      CALL creatg(fidres, "/data/var2d/sv_ky_pj", "singular values of the moment ky/pj cut")
#endif
    ENDIF
    !  var3d group (phi,psi, fluid moments, Ni00, Napjz)
    IF (nsave_3d .GT. 0) THEN
     CALL creatg(fidres, "/data/var3d", "3d profiles")
     CALL creatd(fidres, rank, dims,  "/data/var3d/time",     "Time t*c_s/R")
     CALL creatd(fidres, rank, dims, "/data/var3d/cstep", "iteration number")
     IF (write_phi) CALL creatg(fidres, "/data/var3d/phi", "phi")
     IF (write_phi.AND.EM) CALL creatg(fidres, "/data/var3d/psi", "psi")
    IF (write_Na00) THEN
    CALL creatg(fidres, "/data/var3d/Na00", "gyrocenter density ")
    CALL creatg(fidres, "/data/var3d/Napjz", "pj(z) moment spectrum ")
    ENDIF
    IF (write_dens) THEN
    CALL creatg(fidres, "/data/var3d/dens", "density ")
    ENDIF
    IF (write_fvel) THEN
      CALL creatg(fidres, "/data/var3d/upar", "parallel fluid velocity ")
      CALL creatg(fidres, "/data/var3d/uper", "perpendicular fluid velocity ")
    ENDIF
    IF (write_temp) THEN
      CALL creatg(fidres, "/data/var3d/Tper", "perpendicular temperature ")
      CALL creatg(fidres, "/data/var3d/Tpar", "parallel temperature ")
      CALL creatg(fidres, "/data/var3d/temp", "tiotal temperature ")
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
       CALL creatg(fidres, "/data/var5d/moments", "full moments array")
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
    !                       2.2   2d profiles
    IF (nsave_2d .GT. 0) THEN
      IF (MOD(cstep, nsave_2d) == 0) THEN
        CALL diagnose_2d
      ENDIF
    ENDIF
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
    IF (nsave_5d .GT. 0) THEN
      IF (MOD(cstep, nsave_5d) == 0) THEN
        CALL diagnose_5d
      END IF
    END IF
  !_____________________________________________________________________________
  !                   3.   Final diagnostics
  ELSEIF (kstep .EQ. -1) THEN
    CALL attach(fidres, "/data/input","cpu_time",chrono_runt%ttot)
    ! make a checkpoint at last timestep if not crashed
    IF(.NOT. crashed) THEN
      IF(my_id .EQ. 0) write(*,*) 'Saving last state'
      IF (nsave_5d .GT. 0) CALL diagnose_5d
    ENDIF
    !   Close all diagnostic files
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    CALL closef(fidres)
  END IF
END SUBROUTINE diagnose_full

!! This routine outputs only txt file 0D data (flux and time)
SUBROUTINE diagnose_txtonly(kstep)
  USE basic
  USE diagnostics_par
  USE processing,      ONLY: pflux_x, hflux_x, compute_radial_transport, compute_radial_heatflux
  USE parallel,        ONLY: my_id, comm0
  USE futils,          ONLY: creatf, creatg, creatd, attach, putfile, closef
  IMPLICIT NONE
  INTEGER, INTENT(in) :: kstep
  INTEGER, parameter  :: BUFSIZE = 2
  INTEGER :: rank = 0
  INTEGER :: dims(1) = (/0/)
  IF (kstep .GE. 0) THEN
    ! output the transport in a txt file
    IF ( MOD(cstep, nsave_0d) == 0 ) THEN
      CALL compute_radial_transport
      CALL compute_radial_heatflux
      IF (my_id .EQ. 0) &
        WRITE(1,*) time, pflux_x(1), hflux_x(1)
    END IF
  !! Save the last state
  ELSEIF (kstep .EQ. -1) THEN
    CALL init_outfile(comm0,   resfile0,resfile,fidres)
    CALL creatg(fidres, "/data/var5d", "5d profiles")
    CALL creatd(fidres, rank, dims,  "/data/var5d/time",     "Time t*c_s/R")
    CALL creatd(fidres, rank, dims, "/data/var5d/cstep", "iteration number")
    CALL creatg(fidres, "/data/var5d/moments", "full moments array")
    CALL attach(fidres,"/data/var5d/" , "frames", iframe5d)
    CALL diagnose_5d
    CALL closef(fidres)
  ENDIF
END SUBROUTINE diagnose_txtonly

!!-------------- Auxiliary routines -----------------!!
SUBROUTINE diagnose_0d
  USE basic
  USE futils, ONLY: append, attach, getatt, creatd
  USE diagnostics_par
  USE prec_const
  USE processing
  USE model,   ONLY: Na
  IMPLICIT NONE
  ! Time measurement data
  CALL append(fidres, "/profiler/Tc_rhs",       REAL(chrono_mrhs%ttot,dp),ionode=0)
  CALL append(fidres, "/profiler/Tc_adv_field", REAL(chrono_advf%ttot,dp),ionode=0)
  CALL append(fidres, "/profiler/Tc_clos",      REAL(chrono_clos%ttot,dp),ionode=0)
  CALL append(fidres, "/profiler/Tc_ghost",     REAL(chrono_ghst%ttot,dp),ionode=0)
  CALL append(fidres, "/profiler/Tc_coll",      REAL(chrono_coll%ttot,dp),ionode=0)
  CALL append(fidres, "/profiler/Tc_poisson",   REAL(chrono_pois%ttot,dp),ionode=0)
  CALL append(fidres, "/profiler/Tc_Sapj",      REAL(chrono_sapj%ttot,dp),ionode=0)
  CALL append(fidres, "/profiler/Tc_checkfield",REAL(chrono_chck%ttot,dp),ionode=0)
  CALL append(fidres, "/profiler/Tc_diag",      REAL(chrono_diag%ttot,dp),ionode=0)
  CALL append(fidres, "/profiler/Tc_grad",      REAL(chrono_grad%ttot,dp),ionode=0)
  CALL append(fidres, "/profiler/Tc_nadiab",    REAL(chrono_napj%ttot,dp),ionode=0)
  CALL append(fidres, "/profiler/Tc_step",      REAL(chrono_step%ttot,dp),ionode=0)
#ifdef TEST_SVD
  CALL append(fidres, "/profiler/Tc_CLA",      REAL(chrono_CLA%ttot,dp),ionode=0) 
#endif
  CALL append(fidres, "/profiler/time",                REAL(time,dp),ionode=0)
  ! Processing data
  CALL append(fidres,  "/data/var0d/time",      REAL(time,dp),ionode=0)
  CALL append(fidres, "/data/var0d/cstep", real(cstep,dp),ionode=0)
  CALL getatt(fidres,      "/data/var0d/",       "frames",iframe0d)
  iframe0d=iframe0d+1
  CALL attach(fidres,"/data/var0d/" , "frames", iframe0d)
  ! Ion transport data
  IF (write_gamma) THEN
    CALL compute_radial_transport
    CALL append(fidres, "/data/var0d/gflux_x",REAL(gflux_x,dp),ionode=0)
    CALL append(fidres, "/data/var0d/pflux_x",REAL(pflux_x,dp),ionode=0)
  ENDIF
  IF (write_hf) THEN
    CALL compute_radial_heatflux
    CALL append(fidres, "/data/var0d/hflux_x",REAL(hflux_x,dp),ionode=0)
  ENDIF
END SUBROUTINE diagnose_0d

SUBROUTINE diagnose_2d
  USE prec_const
  USE basic
  USE diagnostics_par
  USE futils, ONLY: putarr, append
#ifdef TEST_SVD
  USE CLA, ONLY: Sf
#endif
  IMPLICIT NONE
  iframe2d=iframe2d+1
  CALL append(fidres,"/data/var2d/time", REAL(time,dp), ionode=0)
  CALL append(fidres,"/data/var2d/cstep",REAL(cstep,dp),ionode=0)
#ifdef TEST_SVD
  WRITE(dset_name, "(A, '/', i6.6)") "/data/var2d/sv_ky_pj/", iframe2d
  CALL putarr(fidres, dset_name, Sf, ionode=0)
#endif
END SUBROUTINE

SUBROUTINE diagnose_3d
  USE basic
  USE futils, ONLY: append, getatt, attach, putarr!, putarrnd
  USE fields, ONLY: phi, psi, moments
  USE array,  ONLY: Napjz,dens,upar,uper,Tpar,Tper,temp
  USE grid, ONLY: CONTAINSp0, ip0,ij0, local_na, total_na,&
                  total_np, total_nj, total_nky, total_nkx, total_nz, &
                  local_np, local_nky, local_nkx, local_nz, &
                  ngz
  USE time_integration, ONLY: updatetlevel
  USE diagnostics_par
  USE prec_const
  USE processing, ONLY: compute_fluid_moments, compute_Napjz_spectrum
  USE model,      ONLY: EM
  USE parallel,   ONLY: manual_3D_bcast
  IMPLICIT NONE
  COMPLEX(xp), DIMENSION(local_na,local_nky,local_nkx,local_nz) :: Na00_
  ! add current time, cstep and frame
  CALL append(fidres,  "/data/var3d/time",           REAL(time,dp),ionode=0)
  CALL append(fidres, "/data/var3d/cstep", real(cstep,dp),ionode=0)
  CALL getatt(fidres,      "/data/var3d/",       "frames",iframe3d)
  iframe3d=iframe3d+1
  CALL attach(fidres,"/data/var3d/" , "frames", iframe3d)
  ! Write current EM fields
  IF (write_phi)        CALL write_field3d_kykxz(phi (:,:,(1+ngz/2):(local_nz+ngz/2)), 'phi')
  IF (write_phi.AND.EM) CALL write_field3d_kykxz(psi (:,:,(1+ngz/2):(local_nz+ngz/2)), 'psi')
  IF (write_Na00) THEN
    CALL compute_Napjz_spectrum
    IF (CONTAINSp0) THEN
      ! gyrocenter density
      Na00_    = moments(:,ip0,ij0,:,:,(1+ngz/2):(local_nz+ngz/2),updatetlevel)
    ELSE
      Na00_    = 0._xp
    ENDIF
    CALL write_field3da_kykxz(Na00_, 'Na00')
    ! <<Napj>x>y spectrum
    CALL write_field3da_pjz(Napjz, 'Napjz')
  ENDIF
  !! Fuid moments
  IF (write_dens .OR. write_fvel .OR. write_temp) &
  CALL compute_fluid_moments
  IF (write_dens) THEN
    CALL write_field3da_kykxz(dens, 'dens')
  ENDIF
  IF (write_fvel) THEN
    CALL write_field3da_kykxz(upar, 'upar')
    CALL write_field3da_kykxz(uper, 'uper')
  ENDIF
  IF (write_temp) THEN
    CALL write_field3da_kykxz(Tpar, 'Tpar')
    CALL write_field3da_kykxz(Tper, 'Tper')
    CALL write_field3da_kykxz(temp, 'temp')
  ENDIF
  CONTAINS
    SUBROUTINE write_field3d_kykxz(field, text)
      USE parallel, ONLY : gather_xyz, num_procs
      IMPLICIT NONE
      COMPLEX(xp), DIMENSION(:,:,:), INTENT(IN) :: field
      CHARACTER(*), INTENT(IN) :: text
      COMPLEX(xp), DIMENSION(total_nky,total_nkx,total_nz) :: field_full
      CHARACTER(50) :: dset_name
      field_full = 0;
      WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var3d", TRIM(text), iframe3d
      IF (num_procs .EQ. 1) THEN ! no data distribution
        CALL putarr(fidres, dset_name, field, ionode=0)
      ELSE
        CALL gather_xyz(field,field_full,local_nky,total_nky,total_nkx,local_nz,total_nz)
        CALL putarr(fidres, dset_name, field_full, ionode=0)
      ENDIF
      END SUBROUTINE write_field3d_kykxz

      SUBROUTINE write_field3da_kykxz(field, text)
        USE parallel, ONLY : gather_xyz, num_procs
        IMPLICIT NONE
        COMPLEX(xp), DIMENSION(:,:,:,:), INTENT(IN) :: field
        CHARACTER(*), INTENT(IN) :: text
        COMPLEX(xp), DIMENSION(total_na,total_nky,total_nkx,total_nz) :: field_full
        COMPLEX(xp), DIMENSION(local_nky,total_nkx,local_nz) :: buff_local
        COMPLEX(xp), DIMENSION(total_nky,total_nkx,total_nz) :: buff_full
        CHARACTER(50) :: dset_name
        INTEGER :: ia
        field_full = 0;
        WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var3d", TRIM(text), iframe3d
        IF (num_procs .EQ. 1) THEN ! no data distribution
          CALL putarr(fidres, dset_name, field, ionode=0)
        ELSE
          DO ia = 1,total_Na
            buff_local = field(ia,:,:,:)
            CALL gather_xyz(buff_local,buff_full,local_nky,total_nky,total_nkx,local_nz,total_nz)
            field_full(ia,:,:,:) = buff_full
          ENDDO
          CALL putarr(fidres, dset_name, field_full, ionode=0)
        ENDIF
        END SUBROUTINE write_field3da_kykxz

    SUBROUTINE write_field3da_pjz(field, text)
      USE parallel, ONLY : gather_pjz, num_procs
      IMPLICIT NONE
      COMPLEX(xp), DIMENSION(:,:,:,:), INTENT(IN) :: field
      CHARACTER(*), INTENT(IN) :: text
      COMPLEX(xp), DIMENSION(total_na,total_np,total_nj,total_nz) :: field_full
      CHARACTER(LEN=50) :: dset_name
      field_full = 0;
      WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var3d", TRIM(text), iframe3d
      IF (num_procs .EQ. 1) THEN ! no data distribution
        CALL putarr(fidres, dset_name, field, ionode=0)
      ELSE
        CALL gather_pjz(field,field_full,total_na,local_np,total_np,total_nj,local_nz,total_nz)
        CALL putarr(fidres, dset_name, field_full, ionode=0)
      ENDIF
      CALL attach(fidres, dset_name, "time", time)
    END SUBROUTINE write_field3da_pjz

END SUBROUTINE diagnose_3d

SUBROUTINE diagnose_5d

  USE basic,  ONLY: time, iframe5d,cstep
  USE futils, ONLY: append, getatt, attach, putarr !,putarrnd
  USE fields, ONLY: moments
  USE grid,   ONLY:total_np, total_nj, total_nky, total_nkx, total_nz, &
                   local_np, local_nj, local_nky, local_nkx, local_nz, &
                   ngp, ngj, ngz, total_na
  USE diagnostics_par
  USE prec_const, ONLY: xp,dp
  IMPLICIT NONE

  CALL append(fidres,  "/data/var5d/time",  REAL(time,dp),ionode=0)
  CALL append(fidres, "/data/var5d/cstep", REAL(cstep,dp),ionode=0)
  CALL getatt(fidres,      "/data/var5d/",       "frames",iframe5d)
  iframe5d=iframe5d+1
  CALL attach(fidres,"/data/var5d/" , "frames", iframe5d)

  IF (write_Napj) THEN
  CALL write_field5d(moments, 'moments')
  ENDIF

  CONTAINS

  SUBROUTINE write_field5d(field, text)
    USE basic,            ONLY: jobnum, dt
    USE futils,           ONLY: attach, putarr!, putarrnd
    USE parallel,         ONLY: gather_pjxyz, num_procs
    USE prec_const,       ONLY: xp
    USE time_integration, ONLY: ntimelevel
    IMPLICIT NONE
    COMPLEX(xp), DIMENSION(total_na,local_np+ngp,local_nj+ngj,local_Nky,local_Nkx,local_nz+ngz,ntimelevel), INTENT(IN) :: field
    CHARACTER(*), INTENT(IN) :: text
    COMPLEX(xp), DIMENSION(total_na,local_np,local_nj,local_nky,local_nkx,local_nz) :: field_sub
    COMPLEX(xp), DIMENSION(total_na,total_np,total_nj,total_nky,total_nkx,total_nz) :: field_full
    CHARACTER(LEN=50) :: dset_name
    field_sub  = field(1:total_na,(1+ngp/2):(local_np+ngp/2),(1+ngj/2):(local_nj+ngj/2),&
                          1:local_nky,1:local_nkx,(1+ngz/2):(local_nz+ngz/2),1)
    field_full = 0;
    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var5d", TRIM(text), iframe5d
    IF (num_procs .EQ. 1) THEN
      CALL putarr(fidres, dset_name, field_sub, ionode=0)
    ELSE!IF(GATHERV_OUTPUT) THEN ! output using one node (gatherv)
      CALL gather_pjxyz(field_sub,field_full,total_na,local_np,total_np,total_nj,local_nky,total_nky,total_nkx,local_nz,total_nz)
      CALL putarr(fidres, dset_name, field_full, ionode=0)
    ! ELSE
    !   CALL putarrnd(fidres, dset_name, field_sub,  (/1,3,5/))
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
  USE prec_const, ONLY: xp
  IMPLICIT NONE
  LOGICAL :: file_exist
  INTEGER :: fid_check, ikx, iky, iz
  CHARACTER(256) :: check_filename
  COMPLEX(xp), DIMENSION(total_nky,total_nkx,total_nz) :: field_to_check
  !! Spit a snapshot of PHI if requested (triggered by creating a file named "check_phi")
  INQUIRE(file='check_phi', exist=file_exist)
  IF( file_exist ) THEN
     IF(my_id.EQ. 0) WRITE(*,*) 'Check file found -> gather phi..'
     CALL gather_xyz(phi(:,:,(1+Ngz/2):(local_nz+Ngz/2)), field_to_check,local_nky,total_nky,total_nkx,local_nz,total_nz)
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
