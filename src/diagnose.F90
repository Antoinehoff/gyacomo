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

  USE prec_const
  IMPLICIT NONE

  INCLUDE 'srcinfo.h'

  INTEGER, INTENT(in) :: kstep
  INTEGER, parameter  :: BUFSIZE = 2
  INTEGER :: rank, dims(1) = (/0/)
  INTEGER :: cp_counter = 0
  CHARACTER(len=256) :: str, fname,test_
  ! putarr(...,pardim=1) does not work for 2D domain decomposition
  ! so we need to gather non 5D data on one proc to output it
  INTEGER     :: parray_e_full(1:pmaxe+1), parray_i_full(1:pmaxi+1)
  INTEGER     :: jarray_e_full(1:jmaxe+1), jarray_i_full(1:jmaxi+1)
  REAL(dp)    :: krarray_full(1:nkr),  kzarray_full(1:nkz)

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
     CALL creatd(fidres, 0, dims, "/profiler/Tc_comm",       "cumulative communication time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_poisson",    "cumulative poisson computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_Sapj",       "cumulative Sapj computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_diag",       "cumulative sym computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_checkfield", "cumulative checkfield computation time")
     CALL creatd(fidres, 0, dims, "/profiler/Tc_step",       "cumulative total step computation time")
     CALL creatd(fidres, 0, dims, "/profiler/time",          "current simulation time")

     ! Build the full grids on process 0 to diagnose it without comm
     IF (my_id .EQ. 0) THEN
       ! P
       DO ip = 1,pmaxe+1; parray_e_full(ip) = (ip-1); END DO
       DO ip = 1,pmaxi+1; parray_i_full(ip) = (ip-1); END DO
       ! J
       DO ij = 1,jmaxe+1; jarray_e_full(ij) = (ij-1); END DO
       DO ij = 1,jmaxi+1; jarray_i_full(ij) = (ij-1); END DO
       ! Kr
       DO ikr = 1,Nkr
         krarray_full(ikr) = REAL(ikr-1,dp) * deltakr
       END DO
       ! Kz
       IF (Nkz .GT. 1) THEN
        DO ikz = 1,Nkz
          kzarray_full(ikz) = deltakz*(MODULO(ikz-1,Nkz/2)-Nkz/2*FLOOR(2.*real(ikz-1)/real(Nkz)))
          if (ikz .EQ. Nz/2+1)     kzarray(ikz) = -kzarray(ikz)
        END DO
      ELSE
        kzarray_full(1) =  0
      endif
     ENDIF

     !  var2d group (electro. pot., Ni00 moment)
     rank = 0
     CALL creatd(fidres, rank, dims,  "/data/var2d/time",     "Time t*c_s/R")
     CALL creatd(fidres, rank, dims, "/data/var2d/cstep", "iteration number")
     CALL creatg(fidres, "/data/grid", "Grid data")
     CALL putarr(fidres, "/data/grid/coordkr", krarray_full(1:nkr),"kr*rho_s0", ionode=0)
     CALL putarr(fidres, "/data/grid/coordkz", kzarray_full(1:nkz),"kz*rho_s0", ionode=0)
     CALL putarr(fidres, "/data/grid/coordp" , parray_e_full(1:pmaxe+1), "p_e", ionode=0)
     CALL putarr(fidres, "/data/grid/coordj" , jarray_e_full(1:jmaxe+1), "j_e", ionode=0)

     IF (nsave_2d .GT. 0) THEN
       CALL creatg(fidres, "/data/var2d/Ne00", "Ne00")
       CALL creatg(fidres, "/data/var2d/Ni00", "Ni00")
       CALL creatg(fidres, "/data/var2d/phi", "phi")
     END IF

     !  var5d group (moments)
     rank = 0
     CALL creatd(fidres, rank, dims,  "/data/var5d/time",     "Time t*c_s/R")
     CALL creatd(fidres, rank, dims, "/data/var5d/cstep", "iteration number")
     IF (nsave_5d .GT. 0) THEN
       CALL creatg(fidres, "/data/var5d/moments_e", "moments_e")
       CALL creatg(fidres, "/data/var5d/moments_i", "moments_i")
      !  CALL creatg(fidres, "/data/var5d/Sepj", "Sepj")
      !  CALL creatg(fidres, "/data/var5d/Sipj", "Sipj")
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
     CALL attach(fidres, TRIM(str),       "Nproc",   num_procs)
     CALL attach(fidres, TRIM(str),       "Np_p" , num_procs_p)
     CALL attach(fidres, TRIM(str),       "Np_kr",num_procs_kr)

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

  !_____________________________________________________________________________
  !                   3.   Final diagnostics

  ELSEIF (kstep .EQ. -1) THEN
     CALL cpu_time(finish)
     CALL attach(fidres, "/data/input","cpu_time",finish-start)

     ! Display computational time cost
     IF (my_id .EQ. 0) CALL display_h_min_s(finish-start)

     !   Close all diagnostic files
     CALL closef(fidres)

  END IF
  
  CALL cpu_time(t1_diag); tc_diag = tc_diag + (t1_diag - t0_diag)

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
  CALL append(fidres, "/profiler/Tc_comm",            tc_comm,ionode=0)
  CALL append(fidres, "/profiler/Tc_step",            tc_step,ionode=0)
  CALL append(fidres, "/profiler/time",                  time,ionode=0)

END SUBROUTINE diagnose_0d


SUBROUTINE diagnose_2d

  USE basic
  USE futils, ONLY: append, getatt, attach, putarrnd
  USE fields
  USE array, ONLY: Ne00, Ni00
  USE grid, ONLY: ikrs,ikre, ikzs,ikze, nkr, nkz, local_nkr, ikr, ikz, ips_e, ips_i
  USE time_integration
  USE diagnostics_par
  USE prec_const
  IMPLICIT NONE

  COMPLEX(dp) :: buffer(ikrs:ikre,ikzs:ikze)
  INTEGER     :: i_, root, world_rank, world_size

  CALL append(fidres,  "/data/var2d/time",           time,ionode=0)
  CALL append(fidres, "/data/var2d/cstep", real(cstep,dp),ionode=0)
  CALL getatt(fidres,      "/data/var2d/",       "frames",iframe2d)
  iframe2d=iframe2d+1
  CALL attach(fidres,"/data/var2d/" , "frames", iframe2d)

  CALL write_field2d(phi (:,:), 'phi')

  IF ( (ips_e .EQ. 1) .AND. (ips_i .EQ. 1) ) THEN
    Ne00(ikrs:ikre,ikzs:ikze) = moments_e(ips_e,1,ikrs:ikre,ikzs:ikze,updatetlevel)
    Ni00(ikrs:ikre,ikzs:ikze) = moments_i(ips_e,1,ikrs:ikre,ikzs:ikze,updatetlevel)
  ENDIF

  root = 0

  !!!!! This is a manual way to do MPI_BCAST !!!!!!!!!!!
  CALL MPI_COMM_RANK(commp,world_rank,ierr)
  CALL MPI_COMM_SIZE(commp,world_size,ierr)

  IF (world_size .GT. 1) THEN
    !! Broadcast phi to the other processes on the same k range (communicator along p)
    IF (world_rank .EQ. root) THEN
      ! Fill the buffer
      DO ikr = ikrs,ikre
        DO ikz = ikzs,ikze
          buffer(ikr,ikz) = Ne00(ikr,ikz)
        ENDDO
      ENDDO
      ! Send it to all the other processes
      DO i_ = 0,num_procs_p-1
        IF (i_ .NE. world_rank) &
        CALL MPI_SEND(buffer, local_nkr * nkz , MPI_DOUBLE_COMPLEX, i_, 0, commp, ierr)
      ENDDO
    ELSE
      ! Recieve buffer from root
      CALL MPI_RECV(buffer, local_nkr * nkz , MPI_DOUBLE_COMPLEX, root, 0, commp, MPI_STATUS_IGNORE, ierr)
      ! Write it in phi
      DO ikr = ikrs,ikre
        DO ikz = ikzs,ikze
          Ne00(ikr,ikz) = buffer(ikr,ikz)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  CALL write_field2d(Ne00(ikrs:ikre,ikzs:ikze), 'Ne00')

    !!!!! This is a manual way to do MPI_BCAST !!!!!!!!!!!
  CALL MPI_COMM_RANK(commp,world_rank,ierr)
  CALL MPI_COMM_SIZE(commp,world_size,ierr)

  IF (world_size .GT. 1) THEN
    !! Broadcast phi to the other processes on the same k range (communicator along p)
    IF (world_rank .EQ. root) THEN
      ! Fill the buffer
      DO ikr = ikrs,ikre
        DO ikz = ikzs,ikze
          buffer(ikr,ikz) = Ni00(ikr,ikz)
        ENDDO
      ENDDO
      ! Send it to all the other processes
      DO i_ = 0,num_procs_p-1
        IF (i_ .NE. world_rank) &
        CALL MPI_SEND(buffer, local_nkr * nkz , MPI_DOUBLE_COMPLEX, i_, 0, commp, ierr)
      ENDDO
    ELSE
      ! Recieve buffer from root
      CALL MPI_RECV(buffer, local_nkr * nkz , MPI_DOUBLE_COMPLEX, root, 0, commp, MPI_STATUS_IGNORE, ierr)
      ! Write it in phi
      DO ikr = ikrs,ikre
        DO ikz = ikzs,ikze
          Ni00(ikr,ikz) = buffer(ikr,ikz)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  CALL write_field2d(Ni00(ikrs:ikre,ikzs:ikze), 'Ni00')

CONTAINS

  SUBROUTINE write_field2d(field, text)
    USE futils, ONLY: attach, putarr
    USE grid, ONLY: ikrs,ikre, ikzs,ikze, nkr, nkz, local_nkr
    USE prec_const
    USE basic, ONLY : commr, num_procs_p, rank_p
    IMPLICIT NONE

    COMPLEX(dp), DIMENSION(ikrs:ikre, ikzs:ikze), INTENT(IN) :: field
    CHARACTER(*), INTENT(IN) :: text
    COMPLEX(dp) :: buffer_dist(ikrs:ikre,ikzs:ikze)
    COMPLEX(dp) :: buffer_full(1:nkr,1:nkz)
    INTEGER     :: scount, rcount
    CHARACTER(LEN=50) :: dset_name

    scount = (ikre-ikrs+1) * (ikze-ikzs+1)
    rcount = scount

    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var2d", TRIM(text), iframe2d
    IF (num_procs .EQ. 1) THEN ! no data distribution
      CALL putarr(fidres, dset_name, field(ikrs:ikre, ikzs:ikze), ionode=0)
    ELSE
      CALL putarrnd(fidres, dset_name, field(ikrs:ikre, ikzs:ikze),  (/1, 1/))
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

   CALL write_field5d_e(moments_e(ips_e:ipe_e,ijs_e:ije_e,:,:,updatetlevel), 'moments_e')
   CALL write_field5d_i(moments_i(ips_i:ipe_i,ijs_i:ije_i,:,:,updatetlevel), 'moments_i')

  !  CALL write_field5d_e(Sepj(ips_e:ipe_e,ijs_e:ije_e,:,:), 'Sepj')
  !  CALL write_field5d_i(Sipj(ips_i:ipe_i,ijs_i:ije_i,:,:), 'Sipj')

 CONTAINS

   SUBROUTINE write_field5d_e(field, text)
     USE futils, ONLY: attach, putarr, putarrnd
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
       CALL putarrnd(fidres, dset_name, field(ips_e:ipe_e,ijs_e:ije_e,ikrs:ikre,ikzs:ikze),  (/1,3/))
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
        CALL putarrnd(fidres, dset_name, field(ips_i:ipe_i,ijs_i:ije_i,ikrs:ikre,ikzs:ikze),  (/1,3/))
      ENDIF
     CALL attach(fidres, dset_name, 'cstep', cstep)
     CALL attach(fidres, dset_name, 'time', time)
     CALL attach(fidres, dset_name, 'jobnum', jobnum)
     CALL attach(fidres, dset_name, 'dt', dt)
     CALL attach(fidres, dset_name, 'iframe2d', iframe2d)
     CALL attach(fidres, dset_name, 'iframe5d', iframe5d)

    END SUBROUTINE write_field5d_i

END SUBROUTINE diagnose_5d
