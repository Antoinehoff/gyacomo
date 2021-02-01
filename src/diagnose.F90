SUBROUTINE diagnose(kstep)
  !   Diagnostics, writing simulation state to disk

  USE basic
  USE grid
  USE diagnostics_par
  USE futils, ONLY: creatf, creatg, creatd, closef, putarr, putfile, attach, openf
  USE model
  USE initial_par
  USE fields
  USE time_integration

  USE prec_const
  IMPLICIT NONE

  INCLUDE 'srcinfo.h'

  INTEGER, INTENT(in) :: kstep
  INTEGER, parameter  :: BUFSIZE = 2
  INTEGER :: rank, dims(1) = (/0/)
  INTEGER :: cp_counter = 0
  CHARACTER(len=256) :: str, fname,test_


  !_____________________________________________________________________________
  !                   1.   Initial diagnostics

  IF ((kstep .EQ. 0)) THEN
    ! Writing output filename
     WRITE(resfile,'(a,a1,i2.2,a3)') TRIM(resfile0),'_',jobnum,'.h5'
     !                      1.1   Initial run
     ! Main output file creation
     IF (write_doubleprecision) THEN
        CALL creatf(resfile, fidres, real_prec='d', mpicomm=MPI_COMM_WORLD)
     ELSE
        CALL creatf(resfile, fidres, mpicomm=MPI_COMM_WORLD)
     END IF
     IF (my_id .EQ. 0) WRITE(*,'(3x,a,a)') TRIM(resfile), ' created'

     ! Checkpoint file creation
     IF (nsave_cp .GT. 0) THEN
       WRITE(rstfile,'(a,a1,i2.2,a3)') TRIM(rstfile0),'_',jobnum,'.h5'
       CALL creatf(rstfile, fidrst, real_prec='d', mpicomm=MPI_COMM_WORLD)
       CALL creatg(fidrst, '/Basic', 'Basic data')
       CALL creatg(fidrst, '/Basic/moments_e', 'electron moments')
       CALL creatg(fidrst, '/Basic/moments_i', 'ion moments')
       CALL creatg(fidrst, '/Basic/phi', 'ES potential')
       ! Attaching informations about moments
       CALL attach(fidrst,"/Basic/moments_e/" , "pmaxe", pmaxe)
       CALL attach(fidrst,"/Basic/moments_e/" , "jmaxe", jmaxe)
       CALL attach(fidrst,"/Basic/moments_e/" , "Trunc", CLOS)
       CALL attach(fidrst,"/Basic/moments_i/" , "pmaxi", pmaxi)
       CALL attach(fidrst,"/Basic/moments_i/" , "jmaxi", jmaxi)
       CALL attach(fidrst,"/Basic/moments_i/" , "Trunc", CLOS)

       IF (my_id .EQ. 0) WRITE(*,'(3x,a,a)') TRIM(rstfile), ' created'
       CALL flush(6)
     ELSE
       IF (my_id .EQ. 0) WRITE(*,'(3x,a,a)') 'No checkpoint'
     ENDIF

     !  Data group
     CALL creatg(fidres, "/data", "data")
     CALL creatg(fidres, "/data/var2d", "2d profiles")
     CALL creatg(fidres, "/data/var5d", "5d profiles")


     ! Initialize counter of number of saves for each category
     IF (cstep==0) THEN
         iframe2d=0
     ENDIF
     CALL attach(fidres,"/data/var2d/" , "frames", iframe2d)
     IF (cstep==0) THEN
         iframe5d=0
     END IF
     CALL attach(fidres,"/data/var5d/" , "frames", iframe5d)

     !  File group
     CALL creatg(fidres, "/files", "files")
     CALL attach(fidres, "/files",  "jobnum", jobnum)

     ! Profiler time measurement
     CALL creatg(fidres, "/profiler", "performance analysis")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_rhs",        "cumulative rhs computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_adv_field",  "cumulative adv. fields computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_poisson",    "cumulative poisson computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_Sapj",       "cumulative Sapj computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_diag",        "cumulative sym computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_checkfield", "cumulative checkfield computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_step",       "cumulative total step computation time")
     CALL creatd(fidres, 0, dims, "/profiler/time",          "current simulation time")

     !  var2d group (electro. pot., Ni00 moment)
     rank = 0
     CALL creatd(fidres, rank, dims,  "/data/var2d/time",     "Time t*c_s/R")
     CALL creatd(fidres, rank, dims, "/data/var2d/cstep", "iteration number")

     IF (nsave_2d .GT. 0) THEN
       CALL creatg(fidres, "/data/var2d/Ne00", "Ne00")
       CALL creatg(fidres, "/data/var2d/Ni00", "Ni00")
       CALL creatg(fidres, "/data/var2d/phi", "phi")
       IF (num_procs .EQ. 1) THEN
         CALL putarr(fidres, "/data/var2d/Ne00/coordkr", krarray(ikrs:ikre), "kr*rho_s0", ionode=0)
         CALL putarr(fidres, "/data/var2d/Ni00/coordkr", krarray(ikrs:ikre), "kr*rho_s0", ionode=0)
         CALL putarr(fidres, "/data/var2d/phi/coordkr", krarray(ikrs:ikre), "kr*rho_s0", ionode=0)
       ELSE
         CALL putarr(fidres, "/data/var2d/Ne00/coordkr", krarray(ikrs:ikre), "kr*rho_s0", pardim=1)
         CALL putarr(fidres, "/data/var2d/Ni00/coordkr", krarray(ikrs:ikre), "kr*rho_s0", pardim=1)
         CALL putarr(fidres, "/data/var2d/phi/coordkr", krarray(ikrs:ikre), "kr*rho_s0", pardim=1)
       ENDIF
       CALL putarr(fidres, "/data/var2d/Ne00/coordkz", kzarray(ikzs:ikze), "kz*rho_s0", ionode=0)
       CALL putarr(fidres, "/data/var2d/Ni00/coordkz", kzarray(ikzs:ikze), "kz*rho_s0", ionode=0)
       CALL putarr(fidres, "/data/var2d/phi/coordkz",  kzarray(ikzs:ikze), "kz*rho_s0", ionode=0)
     END IF

     !  var5d group (moments)
     rank = 0
     CALL creatd(fidres, rank, dims,  "/data/var5d/time",     "Time t*c_s/R")
     CALL creatd(fidres, rank, dims, "/data/var5d/cstep", "iteration number")
     IF (nsave_5d .GT. 0) THEN
       CALL creatg(fidres, "/data/var5d/moments_e", "moments_e")
       CALL creatg(fidres, "/data/var5d/moments_i", "moments_i")
       CALL creatg(fidres, "/data/var5d/Sepj", "Sepj")
       CALL creatg(fidres, "/data/var5d/Sipj", "Sipj")

       CALL putarr(fidres,  "/data/var5d/moments_e/coordp", parray_e(ips_e:ipe_e),       "p_e", ionode=0)
       CALL putarr(fidres,  "/data/var5d/moments_e/coordj", jarray_e(ijs_e:ije_e),       "j_e", ionode=0)
       CALL putarr(fidres,  "/data/var5d/moments_i/coordp", parray_i(ips_i:ipe_i),       "p_i", ionode=0)
       CALL putarr(fidres,  "/data/var5d/moments_i/coordj", jarray_i(ijs_i:ije_i),       "j_i", ionode=0)
       CALL putarr(fidres,  "/data/var5d/Sepj/coordp", parray_e(ips_e:ipe_e),       "p_e", ionode=0)
       CALL putarr(fidres,  "/data/var5d/Sepj/coordj", jarray_e(ijs_e:ije_e),       "j_e", ionode=0)
       CALL putarr(fidres,  "/data/var5d/Sipj/coordp", parray_i(ips_i:ipe_i),       "p_i", ionode=0)
       CALL putarr(fidres,  "/data/var5d/Sipj/coordj", jarray_i(ijs_i:ije_i),       "j_i", ionode=0)
       IF (num_procs .EQ. 1) THEN
         CALL putarr(fidres, "/data/var5d/moments_e/coordkr",    krarray(ikrs:ikre), "kr*rho_s0", ionode=0)
         CALL putarr(fidres, "/data/var5d/moments_i/coordkr",    krarray(ikrs:ikre), "kr*rho_s0", ionode=0)
         CALL putarr(fidres, "/data/var5d/Sepj/coordkr",    krarray(ikrs:ikre), "kr*rho_s0", ionode=0)
         CALL putarr(fidres, "/data/var5d/Sipj/coordkr",    krarray(ikrs:ikre), "kr*rho_s0", ionode=0)
       ELSE
         CALL putarr(fidres, "/data/var5d/moments_e/coordkr",    krarray(ikrs:ikre), "kr*rho_s0", pardim=1)
         CALL putarr(fidres, "/data/var5d/moments_i/coordkr",    krarray(ikrs:ikre), "kr*rho_s0", pardim=1)
         CALL putarr(fidres, "/data/var5d/Sepj/coordkr",    krarray(ikrs:ikre), "kr*rho_s0", pardim=1)
         CALL putarr(fidres, "/data/var5d/Sipj/coordkr",    krarray(ikrs:ikre), "kr*rho_s0", pardim=1)
       ENDIF
       CALL putarr(fidres, "/data/var5d/moments_e/coordkz",    kzarray(ikzs:ikze), "kz*rho_s0", ionode=0)
       CALL putarr(fidres, "/data/var5d/moments_i/coordkz",    kzarray(ikzs:ikze), "kz*rho_s0", ionode=0)
       CALL putarr(fidres, "/data/var5d/Sepj/coordkz",    kzarray(ikzs:ikze), "kz*rho_s0", ionode=0)
       CALL putarr(fidres, "/data/var5d/Sipj/coordkz",    kzarray(ikzs:ikze), "kz*rho_s0", ionode=0)
     END IF

     !  Add input namelist variables as attributes of /data/input, defined in srcinfo.h
     IF (my_id .EQ. 0) WRITE(*,*) 'VERSION=', VERSION
     IF (my_id .EQ. 0) WRITE(*,*)  'BRANCH=', BRANCH
     IF (my_id .EQ. 0) WRITE(*,*)  'AUTHOR=', AUTHOR
     IF (my_id .EQ. 0) WRITE(*,*)    'HOST=', HOST

     WRITE(str,'(a,i2.2)') "/data/input"

     rank=0
     CALL creatd(fidres, rank,dims,TRIM(str),'Input parameters')
     CALL attach(fidres, TRIM(str),     "version",  VERSION) !defined in srcinfo.h
     CALL attach(fidres, TRIM(str),      "branch",   BRANCH) !defined in srcinfo.h
     CALL attach(fidres, TRIM(str),      "author",   AUTHOR) !defined in srcinfo.h
     CALL attach(fidres, TRIM(str),    "execdate", EXECDATE) !defined in srcinfo.h
     CALL attach(fidres, TRIM(str),        "host",     HOST) !defined in srcinfo.h
     CALL attach(fidres, TRIM(str),  "start_time",     time)
     CALL attach(fidres, TRIM(str), "start_cstep",    cstep-1)
     CALL attach(fidres, TRIM(str), "start_iframe2d", iframe2d)
     CALL attach(fidres, TRIM(str), "start_iframe5d", iframe5d)
     CALL attach(fidres, TRIM(str),          "dt",       dt)
     CALL attach(fidres, TRIM(str),        "tmax",     tmax)
     CALL attach(fidres, TRIM(str),        "nrun",     nrun)
     CALL attach(fidres, TRIM(str),    "cpu_time",       -1)
     CALL attach(fidres, TRIM(str),       "Nproc",num_procs)

     CALL grid_outputinputs(fidres, str)

     CALL output_par_outputinputs(fidres, str)

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
     WRITE(*,"(F5.0,A,F5.0)") time,"/",tmax
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
     IF (nsave_2d .GT. 0) THEN
        IF (MOD(cstep, nsave_2d) == 0) THEN
           CALL diagnose_2d
        END IF
     END IF

     !                       2.4   3d profiles
     IF (nsave_5d .GT. 0) THEN
        IF (MOD(cstep, nsave_5d) == 0) THEN
           CALL diagnose_5d
        END IF
     END IF

     !                       2.5   Backups
     IF (nsave_cp .GT. 0) THEN
       IF (MOD(cstep, nsave_cp) == 0) THEN
         CALL checkpoint_save(cp_counter)
         cp_counter = cp_counter + 1
       ENDIF
     ENDIF
  !_____________________________________________________________________________
  !                   3.   Final diagnostics

  ELSEIF (kstep .EQ. -1) THEN
     CALL cpu_time(finish)
     CALL attach(fidres, "/data/input","cpu_time",finish-start)

     ! Display computational time cost
     IF (my_id .EQ. 0) CALL display_h_min_s(finish-start)

     !   Close all diagnostic files
     CALL closef(fidres)
     IF ((nsave_cp .GT. 0) .AND. (.NOT. crashed)) THEN
      CALL checkpoint_save(cp_counter)
      CALL closef(fidrst)
     ENDIF

  END IF

END SUBROUTINE diagnose


SUBROUTINE diagnose_0d

  USE basic
  USE futils, ONLY: append
  USE diagnostics_par
  USE prec_const

  IMPLICIT NONE

  CALL append(fidres, "/profiler/Tc_rhs",              tc_rhs,ionode=0)
  CALL append(fidres, "/profiler/Tc_adv_field",  tc_adv_field,ionode=0)
  CALL append(fidres, "/profiler/Tc_poisson",      tc_poisson,ionode=0)
  CALL append(fidres, "/profiler/Tc_Sapj",            tc_Sapj,ionode=0)
  CALL append(fidres, "/profiler/Tc_diag",            tc_diag,ionode=0)
  CALL append(fidres, "/profiler/Tc_checkfield",tc_checkfield,ionode=0)
  CALL append(fidres, "/profiler/Tc_step",            tc_step,ionode=0)
  CALL append(fidres, "/profiler/time",                  time,ionode=0)

END SUBROUTINE diagnose_0d


SUBROUTINE diagnose_2d

  USE basic
  USE futils, ONLY: append, getatt, attach, putarrnd
  USE fields
  USE time_integration
  USE diagnostics_par
  USE prec_const
  IMPLICIT NONE

  CALL append(fidres,  "/data/var2d/time",           time,ionode=0)
  CALL append(fidres, "/data/var2d/cstep", real(cstep,dp),ionode=0)
  CALL getatt(fidres,      "/data/var2d/",       "frames",iframe2d)
  iframe2d=iframe2d+1
  CALL attach(fidres,"/data/var2d/" , "frames", iframe2d)

  CALL write_field2d(phi(:,:), 'phi')
  CALL write_field2d(moments_e(1,1,:,:,updatetlevel), 'Ne00')
  CALL write_field2d(moments_i(1,1,:,:,updatetlevel), 'Ni00')

CONTAINS

  SUBROUTINE write_field2d(field, text)
    USE futils, ONLY: attach, putarr
    USE grid, ONLY: ikrs,ikre, ikzs,ikze
    USE prec_const
    IMPLICIT NONE

    COMPLEX(dp), DIMENSION(ikrs:ikre, ikzs:ikze), INTENT(IN) :: field
    CHARACTER(*), INTENT(IN) :: text

    CHARACTER(LEN=50) :: dset_name

    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var2d", TRIM(text), iframe2d
    IF (num_procs .EQ. 1) THEN
      CALL putarr(fidres, dset_name, field(ikrs:ikre, ikzs:ikze), ionode=0)
    ELSE
      CALL putarr(fidres, dset_name, field(ikrs:ikre, ikzs:ikze), pardim=1)
    ENDIF

    CALL attach(fidres, dset_name, "time", time)

  END SUBROUTINE write_field2d

END SUBROUTINE diagnose_2d

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

   CALL write_field5d_e(moments_e(:,:,:,:,updatetlevel), 'moments_e')
   CALL write_field5d_i(moments_i(:,:,:,:,updatetlevel), 'moments_i')

   CALL write_field5d_e(Sepj, 'Sepj')
   CALL write_field5d_i(Sipj, 'Sipj')

 CONTAINS

   SUBROUTINE write_field5d_e(field, text)
     USE futils, ONLY: attach, putarr
     USE grid,   ONLY: ips_e,ipe_e, ijs_e,ije_e, ikrs,ikre, ikzs,ikze
     USE prec_const
     IMPLICIT NONE

     COMPLEX(dp), DIMENSION(ips_e:ipe_e,ijs_e:ije_e,ikrs:ikre,ikzs:ikze), INTENT(IN) :: field
     CHARACTER(*), INTENT(IN) :: text

     CHARACTER(LEN=50) :: dset_name

     WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var5d", TRIM(text), iframe5d
     IF (num_procs .EQ. 1) THEN
       CALL putarr(fidres, dset_name, field(ips_e:ipe_e,ijs_e:ije_e,ikrs:ikre,ikzs:ikze), ionode=0)
     ELSE
       CALL putarr(fidres, dset_name, field(ips_e:ipe_e,ijs_e:ije_e,ikrs:ikre,ikzs:ikze), pardim=3)
     ENDIF
     CALL attach(fidres, dset_name, "time", time)

   END SUBROUTINE write_field5d_e

   SUBROUTINE write_field5d_i(field, text)
      USE futils, ONLY: attach, putarr
      USE grid, ONLY: ips_i,ipe_i, ijs_i,ije_i, ikrs,ikre, ikzs,ikze
      USE prec_const
      IMPLICIT NONE

      COMPLEX(dp), DIMENSION(ips_i:ipe_i,ijs_i:ije_i,ikrs:ikre,ikzs:ikze), INTENT(IN) :: field
      CHARACTER(*), INTENT(IN) :: text

      CHARACTER(LEN=50) :: dset_name

      WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var5d", TRIM(text), iframe5d
      IF (num_procs .EQ. 1) THEN
        CALL putarr(fidres, dset_name, field(ips_i:ipe_i,ijs_i:ije_i,ikrs:ikre,ikzs:ikze), ionode=0)
      ELSE
        CALL putarr(fidres, dset_name, field(ips_i:ipe_i,ijs_i:ije_i,ikrs:ikre,ikzs:ikze), pardim=3)
      ENDIF
      CALL attach(fidres, dset_name, "time", time)

    END SUBROUTINE write_field5d_i

END SUBROUTINE diagnose_5d

SUBROUTINE checkpoint_save(cp_step)
  USE basic
  USE grid, ONLY: ips_i,ipe_i, ijs_i,ije_i, ips_e,ipe_e, ijs_e,ije_e, ikrs,ikre, ikzs,ikze
  USE diagnostics_par
  USE futils, ONLY: putarr,attach
  USE model
  USE initial_par
  USE fields
  USE time_integration
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: cp_step
  CHARACTER(LEN=50) :: dset_name

  ! Write state of system to restart file
  WRITE(dset_name, "(A, '/', i6.6)") "/Basic/moments_e", cp_step
  IF (num_procs .EQ. 1) THEN
    CALL putarr(fidrst, dset_name, moments_e(ips_e:ipe_e,ijs_e:ije_e,&
                                                      ikrs:ikre,ikzs:ikze,1), ionode=0)
  ELSE
    CALL putarr(fidrst, dset_name, moments_e(ips_e:ipe_e,ijs_e:ije_e,&
                                                      ikrs:ikre,ikzs:ikze,1), pardim=3)
  ENDIF

  CALL attach(fidrst, dset_name, 'cstep', cstep)
  CALL attach(fidrst, dset_name, 'time', time)
  CALL attach(fidrst, dset_name, 'jobnum', jobnum)
  CALL attach(fidrst, dset_name, 'dt', dt)
  CALL attach(fidrst, dset_name, 'iframe2d', iframe2d)
  CALL attach(fidrst, dset_name, 'iframe5d', iframe5d)

  WRITE(dset_name, "(A, '/', i6.6)") "/Basic/moments_i", cp_step
  IF (num_procs .EQ. 1) THEN
    CALL putarr(fidrst, dset_name, moments_i(ips_i:ipe_i,ijs_i:ije_i,&
                                                      ikrs:ikre,ikzs:ikze,1), ionode=0)
  ELSE
    CALL putarr(fidrst, dset_name, moments_i(ips_i:ipe_i,ijs_i:ije_i,&
                                                      ikrs:ikre,ikzs:ikze,1), pardim=3)
  ENDIF

  CALL attach(fidrst, dset_name, 'cstep', cstep)
  CALL attach(fidrst, dset_name, 'time', time)
  CALL attach(fidrst, dset_name, 'jobnum', jobnum)
  CALL attach(fidrst, dset_name, 'dt', dt)
  CALL attach(fidrst, dset_name, 'iframe2d', iframe2d)
  CALL attach(fidrst, dset_name, 'iframe5d', iframe5d)

  ! Write state of system to restart file
  WRITE(dset_name, "(A, '/', i6.6)") "/Basic/phi", cp_step
  IF (num_procs .EQ. 1) THEN
    CALL putarr(fidrst, dset_name, phi(ikrs:ikre,ikzs:ikze), ionode=0)
  ELSE
    CALL putarr(fidrst, dset_name, phi(ikrs:ikre,ikzs:ikze), pardim=1)
  ENDIF

  CALL attach(fidrst, dset_name, 'cstep', cstep)
  CALL attach(fidrst, dset_name, 'time', time)
  CALL attach(fidrst, dset_name, 'jobnum', jobnum)
  CALL attach(fidrst, dset_name, 'dt', dt)
  CALL attach(fidrst, dset_name, 'iframe2d', iframe2d)
  CALL attach(fidrst, dset_name, 'iframe5d', iframe5d)

  IF (my_id .EQ. 0) THEN
  WRITE(*,'(3x,a)') "Checkpoint file "//TRIM(rstfile)//" updated"
  ENDIF

END SUBROUTINE checkpoint_save
