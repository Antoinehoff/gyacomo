SUBROUTINE diagnose(kstep)
  !   Diagnostics, writing simulation state to disk

  USE basic 
  USE fourier_grid  
  USE diagnostics_par
  USE futils, ONLY: creatf, creatg, creatd, closef, putarr, putfile, attach
  USE model
  USE initial_par
  USE fields
  USE time_integration

  use prec_const
  IMPLICIT NONE
  
  INCLUDE 'srcinfo.h'

  INTEGER, INTENT(in) :: kstep
  ! INTEGER, parameter :: BUFSIZE = 20
  INTEGER :: rank, dims(1) = (/0/)
  CHARACTER(len=256) :: str, fname

  !________________________________________________________________________________
  !
  !                   1.   Initial diagnostics
  IF (kstep .EQ. 0) THEN
     WRITE(*,'(a)') '   Initial diagnostics'
    
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
     !CALL creatg(fidres, "/data/var0d", "0d history arrays")
     !CALL creatg(fidres, "/data/var1d", "1d profiles")
     CALL creatg(fidres, "/data/var2d", "2d profiles")
     CALL creatg(fidres, "/data/var3d", "3d profiles")


     ! Initialize counter of number of saves for each category
     !IF (cstep==0) THEN 
     !   iframe1d=0
     !END IF
     !CALL attach(fidres,"/data/var1d/" , "frames", iframe1d)
     IF (cstep==0) THEN 
        iframe2d=0
     END IF
     CALL attach(fidres,"/data/var2d/" , "frames", iframe2d)
     IF (cstep==0) THEN 
      iframe3d=0
     END IF
     CALL attach(fidres,"/data/var3d/" , "frames", iframe3d)

     !  File group
     CALL creatg(fidres, "/files", "files")
     CALL attach(fidres, "/files",  "jobnum", jobnum)

     !  var2d group
     rank = 0
     CALL creatd(fidres, rank, dims,  "/data/var2d/time",     "Time t*c_s/R")
     CALL creatd(fidres, rank, dims, "/data/var2d/cstep", "iteration number")
     IF (write_phi) THEN
      CALL creatg(fidres, "/data/var2d/phi", "phi")
      CALL putarr(fidres, "/data/var2d/phi/coordkr", krarray(ikrs:ikre), "kr*rho_s0",ionode=0)     
      CALL putarr(fidres, "/data/var2d/phi/coordkz", kzarray(ikzs:ikze), "kz*rho_s0",ionode=0)     
     END IF

     !  var3d group
     rank = 0
     CALL creatd(fidres, rank, dims,  "/data/var3d/time",     "Time t*c_s/R")
     CALL creatd(fidres, rank, dims, "/data/var3d/cstep", "iteration number")
     IF (write_moments) THEN
        CALL creatg(fidres, "/data/var3d/moments", "moments")
        CALL putarr(fidres, "/data/var3d/moments/coordpj", pjarray(ipjs:ipje),"(Jmaxa+1)*p+j+1",ionode=0)
        CALL putarr(fidres, "/data/var3d/moments/coordkr", krarray(ikrs:ikre),      "kr*rho_s0",ionode=0)
        CALL putarr(fidres, "/data/var3d/moments/coordkz", kzarray(ikzs:ikze),      "kz*rho_s0",ionode=0)
     END IF

     !  Add input namelist variables as attributes of /data/input, defined in srcinfo.h
     WRITE(*,*) 'VERSION=', VERSION
     WRITE(*,*)  'BRANCH=', BRANCH
     WRITE(*,*)  'AUTHOR=', AUTHOR
     WRITE(*,*)    'HOST=', HOST

     IF(jobnum .LE. 99) THEN
       WRITE(str,'(a,i2.2)') "/data/input.",jobnum
     ELSE
       WRITE(str,'(a,i3.2)') "/data/input.",jobnum
     END IF
     rank=0
     CALL creatd(fidres, rank,dims,TRIM(str),'Input parameters')
     CALL attach(fidres, TRIM(str),     "version",  VERSION) !defined in srcinfo.h
     CALL attach(fidres, TRIM(str),      "branch",   BRANCH) !defined in srcinfo.h
     CALL attach(fidres, TRIM(str),      "author",   AUTHOR) !defined in srcinfo.h
     CALL attach(fidres, TRIM(str),    "execdate", EXECDATE) !defined in srcinfo.h
     CALL attach(fidres, TRIM(str),        "host",     HOST) !defined in srcinfo.h
     CALL attach(fidres, TRIM(str),  "start_time",     time)
     CALL attach(fidres, TRIM(str), "start_cstep",    cstep)
     CALL attach(fidres, TRIM(str),          "dt",       dt)
     CALL attach(fidres, TRIM(str),        "tmax",     tmax)
     CALL attach(fidres, TRIM(str),        "nrun",     nrun)

     CALL fourier_grid_outputinputs(fidres, str)

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

  END IF


  !________________________________________________________________________________
  !                   2.   Periodic diagnostics
  !
  IF (kstep .GT. 0 .OR. kstep .EQ. 0) THEN

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
        END IF
     END IF

     !                       2.4   3d profiles
     IF (nsave_3d .NE. 0) THEN
        IF (MOD(cstep, nsave_3d) == 0) THEN
           CALL diagnose_3d
        END IF
     END IF

     !________________________________________________________________________________
     !                   3.   Final diagnostics

  ELSEIF (kstep .EQ. -1) THEN
     !   Close all diagnostic files
     CALL closef(fidres)

  END IF

END SUBROUTINE diagnose


SUBROUTINE diagnose_0d

  USE basic
  USE futils, ONLY: append
  USE diagnostics_par
  use prec_const

  IMPLICIT NONE
  WRITE(*,'(a,1x,i7.7,a1,i7.7,20x,a,1pe10.3,10x,a,1pe10.3)') &
          '*** Timestep (this run/total) =', step, '/', cstep, 'Time =', time, 'dt =', dt
  WRITE(*,*)
    
  ! flush stdout of all ranks. Usually only rank 0 should write, but error messages might be written from other ranks as well
  CALL FLUSH(stdout)

  !CALL append(fidres,"/data/var0d/time"                 ,time,           ionode=0)
  !CALL append(fidres,"/data/var0d/cstep"                ,real(cstep,dp), ionode=0)

END SUBROUTINE diagnose_0d


SUBROUTINE diagnose_2d

  USE basic
  USE futils, ONLY: append, getatt, attach, putarrnd
  USE fields
  USE time_integration
  USE diagnostics_par
  use prec_const
  IMPLICIT NONE

  CALL append(fidres,  "/data/var2d/time",           time,ionode=0) 
  CALL append(fidres, "/data/var2d/cstep", real(cstep,dp),ionode=0) 
  CALL getatt(fidres,      "/data/var2d/",       "frames",iframe2d) 
  iframe2d=iframe2d+1
  CALL attach(fidres,"/data/var2d/" , "frames", iframe2d) 

  IF (write_phi) THEN
     CALL write_field2d(phi(:,:), 'phi')
  END IF

CONTAINS

  SUBROUTINE write_field2d(field, text)
    USE futils, ONLY: attach, putarr
    USE fourier_grid, only: ikrs,ikre, ikzs,ikze
    use prec_const
    IMPLICIT NONE

    COMPLEX(dp), DIMENSION(ikrs:ikre, ikzs:ikze), INTENT(IN) :: field
    CHARACTER(*), INTENT(IN) :: text

    CHARACTER(LEN=50) :: dset_name

    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var2d", TRIM(text), iframe2d
    CALL putarr(fidres, dset_name, field(ikrs:ikre, ikzs:ikze),ionode=0)

    CALL attach(fidres, dset_name, "time", time)
    
  END SUBROUTINE write_field2d

END SUBROUTINE diagnose_2d
  
SUBROUTINE diagnose_3d

   USE basic
   USE futils, ONLY: append, getatt, attach, putarrnd
   USE fields
   USE time_integration
   USE diagnostics_par
   use prec_const
   IMPLICIT NONE
 
   CALL append(fidres,  "/data/var3d/time",           time,ionode=0) 
   CALL append(fidres, "/data/var3d/cstep", real(cstep,dp),ionode=0) 
   CALL getatt(fidres,      "/data/var3d/",       "frames",iframe3d) 
   iframe3d=iframe3d+1
   CALL attach(fidres,"/data/var3d/" , "frames", iframe3d) 
 
   IF (write_moments) THEN
      CALL write_field3d(moments(:,:,:,updatetlevel), 'moments')
   END IF
 
 CONTAINS
 
   SUBROUTINE write_field3d(field, text)
     USE futils, ONLY: attach, putarr
     USE fourier_grid, only: ipjs,ipje, ikrs,ikre, ikzs,ikze
     use prec_const
     IMPLICIT NONE
 
     COMPLEX(dp), DIMENSION(ipjs:ipje,ikrs:ikre,ikzs:ikze), INTENT(IN) :: field
     CHARACTER(*), INTENT(IN) :: text
 
     CHARACTER(LEN=50) :: dset_name
 
     WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var3d", TRIM(text), iframe3d
     CALL putarr(fidres, dset_name, field(ipjs:ipje,ikrs:ikre,ikzs:ikze),ionode=0)
 
     CALL attach(fidres, dset_name, "time", time)
     
   END SUBROUTINE write_field3d
 
 END SUBROUTINE diagnose_3d