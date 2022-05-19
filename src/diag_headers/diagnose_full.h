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
INTEGER :: cp_counter = 0
CHARACTER(len=256) :: str,test_
!____________________________________________________________________________
!                   1.   Initial diagnostics

IF ((kstep .EQ. 0)) THEN
  CALL init_outfile(comm0,   resfile0,resfile,fidres)

  ! Profiler time measurement
  CALL creatg(fidres, "/profiler", "performance analysis")
  CALL creatd(fidres, 0, dims, "/profiler/Tc_rhs",        "cumulative rhs computation time")
  CALL creatd(fidres, 0, dims, "/profiler/Tc_adv_field",  "cumulative adv. fields computation time")
  CALL creatd(fidres, 0, dims, "/profiler/Tc_clos",       "cumulative closure computation time")
  CALL creatd(fidres, 0, dims, "/profiler/Tc_ghost",       "cumulative communication time")
  CALL creatd(fidres, 0, dims, "/profiler/Tc_coll",       "cumulative collision computation time")
  CALL creatd(fidres, 0, dims, "/profiler/Tc_poisson",    "cumulative poisson computation time")
  CALL creatd(fidres, 0, dims, "/profiler/Tc_Sapj",       "cumulative Sapj computation time")
  CALL creatd(fidres, 0, dims, "/profiler/Tc_checkfield", "cumulative checkfield computation time")
  CALL creatd(fidres, 0, dims, "/profiler/Tc_diag",       "cumulative sym computation time")
  CALL creatd(fidres, 0, dims, "/profiler/Tc_process",    "cumulative process computation time")
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
  CALL putarrnd(fidres, "/data/metric/gyy",            gyy(izs:ize,0:1), (/1, 1, 1/))
  CALL putarrnd(fidres, "/data/metric/gyz",            gyz(izs:ize,0:1), (/1, 1, 1/))
  CALL putarrnd(fidres, "/data/metric/gzz",            gzz(izs:ize,0:1), (/1, 1, 1/))
  CALL putarrnd(fidres, "/data/metric/hatR",          hatR(izs:ize,0:1), (/1, 1, 1/))
  CALL putarrnd(fidres, "/data/metric/hatZ",          hatZ(izs:ize,0:1), (/1, 1, 1/))
  CALL putarrnd(fidres, "/data/metric/hatB",          hatB(izs:ize,0:1), (/1, 1, 1/))
  CALL putarrnd(fidres, "/data/metric/gradxB",      gradxB(izs:ize,0:1), (/1, 1, 1/))
  CALL putarrnd(fidres, "/data/metric/gradyB",      gradyB(izs:ize,0:1), (/1, 1, 1/))
  CALL putarrnd(fidres, "/data/metric/gradzB",      gradzB(izs:ize,0:1), (/1, 1, 1/))
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
   IF (write_phi) THEN
    CALL creatg(fidres, "/data/var3d/phi_gatherv", "phi_gatherv")
   ENDIF

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
     CALL creatg(fidres, "/data/var5d/moments_e", "Sipj")
     CALL creatg(fidres, "/data/var5d/moments_i", "Sepj")
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


   !                       2.3   2d profiles
   IF (nsave_3d .GT. 0) THEN
      IF (MOD(cstep, nsave_3d) == 0) THEN
        CALL diagnose_3d
      END IF
   END IF

   !                       2.4   3d profiles
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
  ENDIF
  IF (write_hf) THEN
    CALL compute_radial_heatflux
    CALL append(fidres, "/data/var0d/hflux_x",hflux_x,ionode=0)
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
  USE parallel, ONLY : gather_xyz
  USE model, ONLY: KIN_E

  IMPLICIT NONE

  INTEGER        :: i_, root, world_rank, world_size
  CHARACTER(256) :: dset_name
  CALL append(fidres,  "/data/var3d/time",           time,ionode=0)
  CALL append(fidres, "/data/var3d/cstep", real(cstep,dp),ionode=0)
  CALL getatt(fidres,      "/data/var3d/",       "frames",iframe3d)
  iframe3d=iframe3d+1
  CALL attach(fidres,"/data/var3d/" , "frames", iframe3d)

  IF (write_phi) CALL write_field3d_kykxz(phi (ikys:ikye,ikxs:ikxe,izs:ize), 'phi')
  IF (write_phi) THEN
    CALL gather_xyz(phi(ikys:ikye,1:Nkx,izs:ize),phi_full(1:Nky,1:Nkx,1:Nz))
    WRITE(dset_name, "(A, '/', i6.6)") "/data/var3d/phi_gatherv", iframe3d
    CALL putarr(fidres, dset_name, phi_full(1:Nky,1:Nkx,1:Nz), ionode=0)
    CALL attach(fidres, dset_name, "time", time)
  ENDIF

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
    IMPLICIT NONE
    COMPLEX(dp), DIMENSION(ikys:ikye,ikxs:ikxe, izs:ize), INTENT(IN) :: field
    CHARACTER(*), INTENT(IN) :: text
    CHARACTER(256) :: dset_name
    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var3d", TRIM(text), iframe3d
    IF (num_procs .EQ. 1) THEN ! no data distribution
      CALL putarr(fidres, dset_name, field(ikys:ikye,ikxs:ikxe, izs:ize), ionode=0)
    ELSE
      CALL putarrnd(fidres, dset_name, field(ikys:ikye,ikxs:ikxe, izs:ize),  (/1, 1, 3/))
    ENDIF
    CALL attach(fidres, dset_name, "time", time)
  END SUBROUTINE write_field3d_kykxz

  SUBROUTINE write_field3d_pjz_i(field, text)
    IMPLICIT NONE
    REAL(dp), DIMENSION(ips_i:ipe_i,ijs_i:ije_i,izs:ize), INTENT(IN) :: field
    CHARACTER(*), INTENT(IN) :: text
    CHARACTER(LEN=50) :: dset_name
    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var3d", TRIM(text), iframe3d
    IF (num_procs .EQ. 1) THEN ! no data distribution
      CALL putarr(fidres, dset_name, field(ips_i:ipe_i,ijs_i:ije_i,izs:ize), ionode=0)
    ELSE
      CALL putarrnd(fidres, dset_name, field(ips_i:ipe_i,ijs_i:ije_i,izs:ize),  (/1, 0, 3/))
    ENDIF
    CALL attach(fidres, dset_name, "time", time)
  END SUBROUTINE write_field3d_pjz_i

  SUBROUTINE write_field3d_pjz_e(field, text)
    IMPLICIT NONE
    REAL(dp), DIMENSION(ips_e:ipe_e,ijs_e:ije_e,izs:ize), INTENT(IN) :: field
    CHARACTER(*), INTENT(IN) :: text
    CHARACTER(LEN=50) :: dset_name
    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var3d", TRIM(text), iframe3d
    IF (num_procs .EQ. 1) THEN ! no data distribution
      CALL putarr(fidres, dset_name, field(ips_e:ipe_e,ijs_e:ije_e,izs:ize), ionode=0)
    ELSE
      CALL putarrnd(fidres, dset_name, field(ips_e:ipe_e,ijs_e:ije_e,izs:ize),  (/1, 0, 3/))
    ENDIF
    CALL attach(fidres, dset_name, "time", time)
  END SUBROUTINE write_field3d_pjz_e

END SUBROUTINE diagnose_3d

SUBROUTINE diagnose_5d

  USE basic
  USE futils, ONLY: append, getatt, attach, putarrnd
  USE fields
  USE array!, ONLY: Sepj, Sipj
  USE grid, ONLY: ips_e,ipe_e, ips_i, ipe_i, &
                 ijs_e,ije_e, ijs_i, ije_i, &
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
    USE grid,   ONLY: ips_e,ipe_e, ijs_e,ije_e, ikxs,ikxe, ikys,ikye, izs,ize
    USE prec_const
    IMPLICIT NONE

    COMPLEX(dp), DIMENSION(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,ikxs:ikxe,izs:ize), INTENT(IN) :: field
    CHARACTER(*), INTENT(IN) :: text

    CHARACTER(LEN=50) :: dset_name

    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var5d", TRIM(text), iframe5d
    IF (num_procs .EQ. 1) THEN
     CALL putarr(fidres, dset_name, field(ips_e:ipe_e,ijs_e:ije_e,ikys:ikye,ikxs:ikxe,izs:ize), ionode=0)
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
    USE grid, ONLY: ips_i,ipe_i, ijs_i,ije_i, ikxs,ikxe, ikys,ikye, izs,ize
    USE prec_const
    IMPLICIT NONE

    COMPLEX(dp), DIMENSION(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,ikxs:ikxe,izs:ize), INTENT(IN) :: field
    CHARACTER(*), INTENT(IN) :: text

    CHARACTER(LEN=50) :: dset_name

    WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var5d", TRIM(text), iframe5d
    IF (num_procs .EQ. 1) THEN
      CALL putarr(fidres, dset_name, field(ips_i:ipe_i,ijs_i:ije_i,ikys:ikye,ikxs:ikxe,izs:ize), ionode=0)
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
