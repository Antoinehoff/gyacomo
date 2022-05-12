SUBROUTINE diagnose_fields(kstep)
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
  INTEGER :: rank = 0
  INTEGER :: dims(1) = (/0/)
  !____________________________________________________________________________!
  IF ((kstep .EQ. 0)) THEN
   !  fields result file (electro. pot., Ni00 moment, fluid moments etc.)
   IF (nsave_3d .GT. 0) THEN
    CALL init_outfile(comm_kyz,fldfile0,fldfile,fidfld)
    CALL creatd(fidfld, rank, dims,  "/data/time",     "Time t*c_s/R")
    CALL creatd(fidfld, rank, dims, "/data/cstep", "iteration number")

    IF (write_phi) CALL creatg(fidfld, "/data/phi", "phi")

    IF (write_Na00) THEN
     IF(KIN_E)&
     CALL creatg(fidfld, "/data/Ne00", "Ne00")
     CALL creatg(fidfld, "/data/Ni00", "Ni00")
    ENDIF

    IF (write_dens) THEN
      IF(KIN_E)&
     CALL creatg(fidfld, "/data/dens_e", "dens_e")
     CALL creatg(fidfld, "/data/dens_i", "dens_i")
    ENDIF

    IF (write_fvel) THEN
      IF(KIN_E) THEN
      CALL creatg(fidfld, "/data/upar_e", "upar_e")
      CALL creatg(fidfld, "/data/uper_e", "uper_e")
      ENDIF
      CALL creatg(fidfld, "/data/upar_i", "upar_i")
      CALL creatg(fidfld, "/data/uper_i", "uper_i")
    ENDIF

    IF (write_temp) THEN
      IF(KIN_E) THEN
      CALL creatg(fidfld, "/data/Tper_e", "Tper_e")
      CALL creatg(fidfld, "/data/Tpar_e", "Tpar_e")
      CALL creatg(fidfld, "/data/temp_e", "temp_e")
      ENDIF
      CALL creatg(fidfld, "/data/Tper_i", "Tper_i")
      CALL creatg(fidfld, "/data/Tpar_i", "Tpar_i")
      CALL creatg(fidfld, "/data/temp_i", "temp_i")
    ENDIF

    IF (cstep==0) THEN
      iframe3d=0
    ENDIF
    CALL attach(fidfld,"/data/" , "frames", iframe3d)
   END IF
  ENDIF

  !_____________________________________________________________________________
  !                   2.   Periodic diagnostics
  !
  IF (kstep .GE. 0) THEN

  !                       2.2   3d profiles
  IF (nsave_3d .GT. 0) THEN
     IF (MOD(cstep, nsave_3d) == 0) THEN

       CALL append(fidfld,  "/data/time",           time,ionode=0)
       CALL append(fidfld, "/data/cstep", real(cstep,dp),ionode=0)
       CALL getatt(fidfld,      "/data/",       "frames",iframe3d)
       iframe3d=iframe3d+1
       CALL attach(fidfld,"/data/" , "frames", iframe3d)

       IF (write_phi) CALL write_field3d_kykxz(phi (ikys:ikye,ikxs:ikxe,izs:ize), 'phi')

       IF (write_Na00) THEN
         IF(KIN_E)THEN
         IF (ips_e .EQ. 1) THEN
           Ne00(ikys:ikye,ikxs:ikxe,izs:ize) = moments_e(ips_e,1,ikys:ikye,ikxs:ikxe,izs:ize,updatetlevel)
         ENDIF
         CALL write_field3d_kykxz(Ne00(ikys:ikye,ikxs:ikxe,izs:ize), 'Ne00')
         ENDIF
         IF (ips_i .EQ. 1) THEN
           Ni00(ikys:ikye,ikxs:ikxe,izs:ize) = moments_i(ips_e,1,ikys:ikye,ikxs:ikxe,izs:ize,updatetlevel)
         ENDIF
         CALL write_field3d_kykxz(Ni00(ikys:ikye,ikxs:ikxe,izs:ize), 'Ni00')
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

     END IF
  END IF
  !_____________________________________________________________________________
  !                   3.   Final diagnostics
  ELSEIF (kstep .EQ. -1) THEN
    !   Close diagnostic files
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    CALL closef(fidfld)

  END IF

  CONTAINS
    SUBROUTINE write_field3d_kykxz(field, text)
      IMPLICIT NONE
      COMPLEX(dp), DIMENSION(ikys:ikye,ikxs:ikxe, izs:ize), INTENT(IN) :: field
      CHARACTER(*), INTENT(IN) :: text
      CHARACTER(LEN=50) :: dset_name
      WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data", TRIM(text), iframe3d
      IF (num_procs .EQ. 1) THEN ! no data distribution
        CALL putarr(fidfld, dset_name, field(ikys:ikye,ikxs:ikxe, izs:ize), ionode=0)
      ELSE
        CALL putarrnd(fidfld, dset_name, field(ikys:ikye,ikxs:ikxe, izs:ize),  (/1, 3/))
      ENDIF
      CALL attach(fidfld, dset_name, "time", time)
    END SUBROUTINE write_field3d_kykxz

END SUBROUTINE diagnose_fields
