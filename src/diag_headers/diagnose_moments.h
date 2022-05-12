SUBROUTINE diagnose_moments(kstep)
  USE basic
  USE grid
  USE diagnostics_par
  USE futils,         ONLY: creatg,creatd,append,putarr,putarrnd,attach,closef,getatt
  USE model,          ONLY: KIN_E
  USE array,          ONLY: Sepj,Sipj
  USE fields,         ONLY: moments_e, moments_i
  USE time_integration
  USE utility
  USE prec_const
  IMPLICIT NONE
  INTEGER, INTENT(in) :: kstep
  INTEGER, parameter  :: BUFSIZE = 2
  INTEGER :: rank = 0
  INTEGER :: dims(1) = (/0/)
  !____________________________________________________________________________!
  IF ((kstep .EQ. 0)) THEN

  !  var5d group (moments)
  IF (nsave_5d .GT. 0) THEN
    CALL init_outfile(comm0,   momfile0,momfile,fidmom)

    CALL creatd(fidmom, rank, dims,  "/data/time",     "Time t*c_s/R")
    CALL creatd(fidmom, rank, dims, "/data/cstep", "iteration number")

    IF (write_Napj) THEN
      IF(KIN_E)&
     CALL creatg(fidmom, "/data/moments_e", "moments_e")
     CALL creatg(fidmom, "/data/moments_i", "moments_i")
    ENDIF

    IF (write_Sapj) THEN
      IF(KIN_E)&
     CALL creatg(fidmom, "/data/moments_e", "Sipj")
     CALL creatg(fidmom, "/data/moments_i", "Sepj")
    ENDIF

    IF (cstep==0) THEN
     iframe5d=0
    END IF
    CALL attach(fidmom,"/data/" , "frames", iframe5d)
  END IF

  ENDIF

  !_____________________________________________________________________________
  !                   2.   Periodic diagnostics
  !
  IF ((kstep .GE. 0) .OR. ((kstep .EQ. -1) .AND. (.NOT. CRASHED))) THEN

  IF((kstep .EQ. -1) .AND. (.NOT. CRASHED .AND. (my_id .EQ. 0))) &
    write(*,*) 'Saving last state'

  !                       2.3   5d profiles
    IF (MOD(cstep, nsave_5d) == 0) THEN
     CALL append(fidmom,  "/data/time",           time,ionode=0)
     CALL append(fidmom, "/data/cstep", real(cstep,dp),ionode=0)
     CALL getatt(fidmom,      "/data/",       "frames",iframe5d)
     iframe5d=iframe5d+1
     CALL attach(fidmom,"/data/" , "frames", iframe5d)

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
    END IF
  END IF
  !_____________________________________________________________________________
  !                   3.   Final diagnostics
  IF (kstep .EQ. -1) THEN
    !   Close diagnostic files
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    CALL closef(fidmom)
  END IF

  CONTAINS
    SUBROUTINE write_field5d_e(field, text)
      USE futils, ONLY: attach, putarr, putarrnd
      USE grid,   ONLY: ips_e,ipe_e, ijs_e,ije_e, ikxs,ikxe, ikys,ikye, izs,ize
      USE prec_const
      IMPLICIT NONE

      COMPLEX(dp), DIMENSION(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,ikxs:ikxe,izs:ize), INTENT(IN) :: field
      CHARACTER(*), INTENT(IN) :: text

      CHARACTER(LEN=50) :: dset_name

      WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data", TRIM(text), iframe5d
      IF (num_procs .EQ. 1) THEN
       CALL putarr(fidmom, dset_name, field(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,ikxs:ikxe,izs:ize), ionode=0)
      ELSE
       CALL putarrnd(fidmom, dset_name, field(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,ikxs:ikxe,izs:ize),  (/1,3,5/))
      ENDIF
      CALL attach(fidmom, dset_name, 'cstep', cstep)
      CALL attach(fidmom, dset_name, 'time', time)
      CALL attach(fidmom, dset_name, 'jobnum', jobnum)
      CALL attach(fidmom, dset_name, 'dt', dt)
      CALL attach(fidmom, dset_name, 'iframe2d', iframe2d)
      CALL attach(fidmom, dset_name, 'iframe5d', iframe5d)

    END SUBROUTINE write_field5d_e

    SUBROUTINE write_field5d_i(field, text)
      USE futils, ONLY: attach, putarr, putarrnd
      USE grid, ONLY: ips_i,ipe_i, ijs_i,ije_i, ikxs,ikxe, ikys,ikye, izs,ize
      USE prec_const
      IMPLICIT NONE

      COMPLEX(dp), DIMENSION(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,ikxs:ikxe,izs:ize), INTENT(IN) :: field
      CHARACTER(*), INTENT(IN) :: text

      CHARACTER(LEN=50) :: dset_name

      WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data", TRIM(text), iframe5d
      IF (num_procs .EQ. 1) THEN
        CALL putarr(fidmom, dset_name, field(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,ikxs:ikxe,izs:ize), ionode=0)
      ELSE
        CALL putarrnd(fidmom, dset_name, field(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,ikxs:ikxe,izs:ize),  (/1,3,5/))
      ENDIF
      CALL attach(fidmom, dset_name, 'cstep', cstep)
      CALL attach(fidmom, dset_name, 'time', time)
      CALL attach(fidmom, dset_name, 'jobnum', jobnum)
      CALL attach(fidmom, dset_name, 'dt', dt)
      CALL attach(fidmom, dset_name, 'iframe2d', iframe2d)
      CALL attach(fidmom, dset_name, 'iframe5d', iframe5d)

    END SUBROUTINE write_field5d_i
END SUBROUTINE diagnose_moments
