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
  CHARACTER(len=256) :: str, fname

  !_____________________________________________________________________________
  !                   1.   Initial diagnostics

  IF ((kstep .EQ. 0)) THEN
     ! jobnum is either 0 from initialization or some integer values read from chkrst(0)
     IF(jobnum .LE. 99) THEN
         WRITE(resfile,'(a,a1,i2.2,a3)') TRIM(resfile0),'_',jobnum,'.h5'
     ELSE
         WRITE(resfile,'(a,a1,i3.2,a3)') TRIM(resfile0),'_',jobnum,'.h5'
     END IF

     !                       1.1   Initial run
     IF (write_doubleprecision) THEN
        CALL creatf(resfile, fidres, real_prec='d')
     ELSE
        CALL creatf(resfile, fidres)
     END IF


     WRITE(*,'(3x,a,a)') TRIM(resfile), ' created'
     call flush(6)


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

     !  var2d group (electro. pot., Ni00 moment)
     rank = 0
     CALL creatd(fidres, rank, dims,  "/data/var2d/time",     "Time t*c_s/R")
     CALL creatd(fidres, rank, dims, "/data/var2d/cstep", "iteration number")

     IF (write_Ni00) THEN
       CALL creatg(fidres, "/data/var2d/Ni00", "Ni00")
       CALL putarr(fidres, "/data/var2d/Ni00/coordkr", krarray(ikrs:ikre), "kr*rho_s0",ionode=0)
       CALL putarr(fidres, "/data/var2d/Ni00/coordkz", kzarray(ikzs:ikze), "kz*rho_s0",ionode=0)
     END IF

     IF (write_phi) THEN
       CALL creatg(fidres, "/data/var2d/phi", "phi")
       CALL putarr(fidres, "/data/var2d/phi/coordkr", krarray(ikrs:ikre), "kr*rho_s0",ionode=0)
       CALL putarr(fidres, "/data/var2d/phi/coordkz", kzarray(ikzs:ikze), "kz*rho_s0",ionode=0)
     END IF

     !  var5d group (moments)
     rank = 0
     CALL creatd(fidres, rank, dims,  "/data/var5d/time",     "Time t*c_s/R")
     CALL creatd(fidres, rank, dims, "/data/var5d/cstep", "iteration number")
     IF (write_moments) THEN
       CALL creatg(fidres, "/data/var5d/moments_e", "moments_e")
       CALL putarr(fidres,  "/data/var5d/moments_e/coordp", parray_e(ips_e:ipe_e),       "p_e",ionode=0)
       CALL putarr(fidres,  "/data/var5d/moments_e/coordj", jarray_e(ijs_e:ije_e),       "j_e",ionode=0)
       CALL putarr(fidres, "/data/var5d/moments_e/coordkr",    krarray(ikrs:ikre), "kr*rho_s0",ionode=0)
       CALL putarr(fidres, "/data/var5d/moments_e/coordkz",    kzarray(ikzs:ikze), "kz*rho_s0",ionode=0)

       CALL creatg(fidres, "/data/var5d/moments_i", "moments_i")
       CALL putarr(fidres,  "/data/var5d/moments_i/coordp", parray_i(ips_i:ipe_i),       "p_i",ionode=0)
       CALL putarr(fidres,  "/data/var5d/moments_i/coordj", jarray_i(ijs_i:ije_i),       "j_i",ionode=0)
       CALL putarr(fidres, "/data/var5d/moments_i/coordkr",    krarray(ikrs:ikre), "kr*rho_s0",ionode=0)
       CALL putarr(fidres, "/data/var5d/moments_i/coordkz",    kzarray(ikzs:ikze), "kz*rho_s0",ionode=0)
     END IF

     IF (write_non_lin) THEN
       CALL creatg(fidres, "/data/var5d/Sepj", "Sepj")
       CALL putarr(fidres,  "/data/var5d/Sepj/coordp", parray_e(ips_e:ipe_e),       "p_e",ionode=0)
       CALL putarr(fidres,  "/data/var5d/Sepj/coordj", jarray_e(ijs_e:ije_e),       "j_e",ionode=0)
       CALL putarr(fidres, "/data/var5d/Sepj/coordkr",    krarray(ikrs:ikre), "kr*rho_s0",ionode=0)
       CALL putarr(fidres, "/data/var5d/Sepj/coordkz",    kzarray(ikzs:ikze), "kz*rho_s0",ionode=0)

       CALL creatg(fidres, "/data/var5d/Sipj", "Sipj")
       CALL putarr(fidres,  "/data/var5d/Sipj/coordp", parray_i(ips_i:ipe_i),       "p_i",ionode=0)
       CALL putarr(fidres,  "/data/var5d/Sipj/coordj", jarray_i(ijs_i:ije_i),       "j_i",ionode=0)
       CALL putarr(fidres, "/data/var5d/Sipj/coordkr",    krarray(ikrs:ikre), "kr*rho_s0",ionode=0)
       CALL putarr(fidres, "/data/var5d/Sipj/coordkz",    kzarray(ikzs:ikze), "kz*rho_s0",ionode=0)
     END IF

     !  Add input namelist variables as attributes of /data/input, defined in srcinfo.h
     WRITE(*,*) 'VERSION=', VERSION
     WRITE(*,*)  'BRANCH=', BRANCH
     WRITE(*,*)  'AUTHOR=', AUTHOR
     WRITE(*,*)    'HOST=', HOST

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
     CALL attach(fidres, TRIM(str),          "dt",       dt)
     CALL attach(fidres, TRIM(str),        "tmax",     tmax)
     CALL attach(fidres, TRIM(str),        "nrun",     nrun)

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

     !                       2.1   0d history arrays
     IF (nsave_0d .NE. 0) THEN
        IF ( MOD(cstep, nsave_0d) == 0 ) THEN
           CALL diagnose_0d
        END IF
     END IF

     !                       2.2   1d profiles
     ! empty in our case

     !                       2.3   2d profiles
     IF (nsave_2d .NE. 0) THEN
        IF (MOD(cstep, nsave_2d) == 0) THEN
           CALL diagnose_2d
           IF (MOD(cstep, nsave_2d*10) == 0) THEN
            WRITE(*,"(F4.0,A,F4.0)") time,"/",tmax
           ENDIF
        END IF
     END IF

     !                       2.4   3d profiles
     IF (nsave_5d .NE. 0) THEN
        IF (MOD(cstep, nsave_5d) == 0) THEN
           CALL diagnose_5d
        END IF
     END IF

     !                       2.5   Backups
     IF (MOD(cstep, nsave_cp) == 0) THEN
       CALL checkpoint_save
     ENDIF

  !_____________________________________________________________________________
  !                   3.   Final diagnostics

  ELSEIF (kstep .EQ. -1) THEN
     CALL cpu_time(finish)
     !   Close all diagnostic files
     CALL closef(fidres)
  END IF

END SUBROUTINE diagnose


SUBROUTINE diagnose_0d

  USE basic
  USE futils, ONLY: append
  USE diagnostics_par
  USE prec_const

  IMPLICIT NONE
  WRITE(*,'(a,1x,i7.7,a1,i7.7,20x,a,1pe10.3,10x,a,1pe10.3)') &
          '*** Timestep (this run/total) =', step, '/', cstep, 'Time =', time, 'dt =', dt
  WRITE(*,*)

  ! flush stdout of all ranks. Usually ONLY rank 0 should write, but error messages might be written from other ranks as well
  CALL FLUSH(stdout)

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

  IF (write_phi) THEN
     CALL write_field2d(phi(:,:), 'phi')
  END IF

  IF (write_Ni00) THEN
     CALL write_field2d(moments_i(1,1,:,:,updatetlevel), 'Ni00')
  END IF
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
    CALL putarr(fidres, dset_name, field(ikrs:ikre, ikzs:ikze),ionode=0)

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

   IF (write_moments) THEN
      CALL write_field5d_e(moments_e(:,:,:,:,updatetlevel), 'moments_e')
      CALL write_field5d_i(moments_i(:,:,:,:,updatetlevel), 'moments_i')
   END IF

   IF (write_non_lin) THEN
      CALL write_field5d_e(Sepj, 'Sepj')
      CALL write_field5d_i(Sipj, 'Sipj')
   END IF

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
     CALL putarr(fidres, dset_name, field(ips_e:ipe_e,ijs_e:ije_e,ikrs:ikre,ikzs:ikze),ionode=0)

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
      CALL putarr(fidres, dset_name, field(ips_i:ipe_i,ijs_i:ije_i,ikrs:ikre,ikzs:ikze),ionode=0)

      CALL attach(fidres, dset_name, "time", time)

    END SUBROUTINE write_field5d_i

END SUBROUTINE diagnose_5d

SUBROUTINE checkpoint_save
  USE basic
  USE grid
  USE diagnostics_par
  USE futils, ONLY: creatf, creatg, creatd, closef, putarr, putfile, attach, openf
  USE model
  USE initial_par
  USE fields
  USE time_integration
  IMPLICIT NONE

  WRITE(rstfile,'(a,a3)') TRIM(rstfile0),'.h5'

  CALL creatf(rstfile, fidrst, real_prec='d')
  CALL creatg(fidrst, '/Basic', 'Basic data')
  CALL attach(fidrst, '/Basic', 'cstep', cstep)
  CALL attach(fidrst, '/Basic', 'time', time)
  CALL attach(fidrst, '/Basic', 'jobnum', jobnum)
  CALL attach(fidrst, '/Basic', 'dt', dt)
  CALL attach(fidrst, '/Basic', 'iframe2d', iframe2d)
  CALL attach(fidrst, '/Basic', 'iframe5d', iframe5d)

  ! Write state of system to restart file
  CALL putarr(fidrst, '/Basic/moments_e', moments_e(ips_e:ipe_e,ijs_e:ije_e,&
                                                    ikrs:ikre,ikzs:ikze,1),ionode=0)
  CALL putarr(fidrst, '/Basic/moments_i', moments_i(ips_i:ipe_i,ijs_i:ije_i,&
                                                    ikrs:ikre,ikzs:ikze,1),ionode=0)
  CALL closef(fidrst)
  WRITE(*,'(3x,a)') "Checkpoint file "//TRIM(rstfile)//" saved!"

END SUBROUTINE checkpoint_save
