SUBROUTINE diagnose_momspectrum(kstep)
  USE basic
  USE grid
  USE diagnostics_par
  USE futils
  USE array
  USE model
  USE initial_par
  USE fields
  USE time_integration
  USE utility
  USE prec_const
  USE collision, ONLY: coll_outputinputs
  USE geometry
  USE processing
  IMPLICIT NONE
  INTEGER, INTENT(in) :: kstep
  INTEGER, parameter  :: BUFSIZE = 2
  INTEGER :: rank = 0, frame
  INTEGER :: dims(1) = (/0/)
  !____________________________________________________________________________!
  IF ((kstep .EQ. 0)) THEN
    !  var3d group (pjz, spectrum of moments)
    IF (nsave_3d .GT. 0) THEN
     CALL init_outfile(comm_pz, mspfile0,mspfile,fidmsp)
     CALL creatd(fidmsp, rank, dims,  "/data/time",     "Time t*c_s/R")
     CALL creatd(fidmsp, rank, dims, "/data/cstep", "iteration number")

     IF (write_Na00) THEN
      IF(KIN_E)&
      CALL creatg(fidmsp, "/data/Nepjz", "Nepjz")
      CALL creatg(fidmsp, "/data/Nipjz", "Nipjz")
     ENDIF

     IF (cstep==0) THEN
       iframe3d=0
     ENDIF
     CALL attach(fidmsp,"/data/" , "frames", iframe3d)
    END IF

  ENDIF

  !_____________________________________________________________________________
  !                   2.   Periodic diagnostics
  !
  IF (kstep .GE. 0) THEN
  !                       2.2   3d profiles
    IF (nsave_3d .GT. 0) THEN
       IF (MOD(cstep, nsave_3d) == 0) THEN
         CALL append(fidmsp,  "/data/time",           time,ionode=0)
         CALL append(fidmsp, "/data/cstep", real(cstep,dp),ionode=0)
         CALL getatt(fidmsp,      "/data/",       "frames",frame)
         frame=frame+1
         CALL attach(fidmsp,"/data/" , "frames", frame)

         IF (write_Na00) THEN
           CALL compute_Napjz_spectrum
           IF(KIN_E) &
           CALL write_field3d_pjz_e(Nepjz(ips_e:ipe_e,ijs_e:ije_e,izs:ize), 'Nepjz', frame)
           CALL write_field3d_pjz_i(Nipjz(ips_i:ipe_i,ijs_i:ije_i,izs:ize), 'Nipjz', frame)
         ENDIF
       END IF
    END IF

  !_____________________________________________________________________________
  !                   3.   Final diagnostics
  ELSEIF (kstep .EQ. -1) THEN
      !   Close diagnostic files
      CALL mpi_barrier(MPI_COMM_WORLD, ierr)
      CALL closef(fidmsp)

  END IF

    CONTAINS
      SUBROUTINE write_field3d_pjz_i(field, text, frame)
        IMPLICIT NONE
        REAL(dp), DIMENSION(ips_i:ipe_i,ijs_i:ije_i,izs:ize), INTENT(IN) :: field
        CHARACTER(*), INTENT(IN) :: text
        INTEGER, INTENT(IN) :: frame
        CHARACTER(LEN=50) :: dset_name
        WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data", TRIM(text), frame
        IF (num_procs .EQ. 1) THEN ! no data distribution
          CALL putarr(fidmsp, dset_name, field(ips_i:ipe_i,ijs_i:ije_i,izs:ize), ionode=0)
        ELSE
          CALL putarrnd(fidmsp, dset_name, field(ips_i:ipe_i,ijs_i:ije_i,izs:ize),  (/1, 3/))
        ENDIF
        CALL attach(fidmsp, dset_name, "time", time)
      END SUBROUTINE write_field3d_pjz_i

      SUBROUTINE write_field3d_pjz_e(field, text, frame)
        IMPLICIT NONE
        REAL(dp), DIMENSION(ips_e:ipe_e,ijs_e:ije_e,izs:ize), INTENT(IN) :: field
        CHARACTER(*), INTENT(IN) :: text
        INTEGER, INTENT(IN) :: frame
        CHARACTER(LEN=50) :: dset_name
        WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data", TRIM(text), frame
        IF (num_procs .EQ. 1) THEN ! no data distribution
          CALL putarr  (fidmsp, dset_name, field(ips_e:ipe_e,ijs_e:ije_e,izs:ize), ionode=0)
        ELSE
          CALL putarrnd(fidmsp, dset_name, field(ips_e:ipe_e,ijs_e:ije_e,izs:ize),  (/1, 3/))
        ENDIF
        CALL attach(fidmsp, dset_name, "time", time)
      END SUBROUTINE write_field3d_pjz_e

END SUBROUTINE diagnose_momspectrum
