MODULE diagnostics
  !   Module for diagnostic parameters

  USE prec_const
  USE h5fortran
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PUBLIC, PROTECTED :: write_doubleprecision = .TRUE.
  LOGICAL, PUBLIC, PROTECTED :: write_gamma = .TRUE.  ! output particle transport and heat flux
  LOGICAL, PUBLIC, PROTECTED :: write_hf    = .TRUE.   ! output particle transport and heat flux
  LOGICAL, PUBLIC, PROTECTED :: write_phi   = .TRUE. 
  LOGICAL, PUBLIC, PROTECTED :: write_Na00  = .TRUE. 
  LOGICAL, PUBLIC, PROTECTED :: write_Napj  = .TRUE. 
  LOGICAL, PUBLIC, PROTECTED :: write_Sapj  = .TRUE. 
  LOGICAL, PUBLIC, PROTECTED :: write_dens  = .TRUE. 
  LOGICAL, PUBLIC, PROTECTED :: write_fvel  = .TRUE. 
  LOGICAL, PUBLIC, PROTECTED :: write_temp  = .TRUE. 

  INTEGER, PUBLIC, PROTECTED :: nsave_0d, nsave_1d, nsave_2d, nsave_3d, nsave_5d ! save data every n step
  REAL,    PUBLIC, PROTECTED :: dtsave_0d, dtsave_1d, dtsave_2d, dtsave_3d, dtsave_5d ! save data every dt time unit
  ! Change diagnostic mode (full/txtonly)
  CHARACTER(len=256), PUBLIC, PROTECTED :: diag_mode = "full"
  !  HDF5 file
  CHARACTER(len=256), PUBLIC :: resfile,resfile0 = "outputs"            ! Head of main result file name
  CHARACTER(len=256), PUBLIC :: momfile,momfile0 = "moments"   ! Head of the moment spectrum file (N_a(p,j,z))
  CHARACTER(len=256), PUBLIC :: mspfile,mspfile0 = "moments_spectrum"   ! Head of the moment spectrum file (N_a(p,j,z))
  CHARACTER(len=256), PUBLIC :: fldfile,fldfile0 = "fields"             ! Head of field (phi,A)
  CHARACTER(len=256), PUBLIC :: ttrfile,ttrfile0 = "time_traces"        ! Head of time traces (gamma_x,Q_x)
  CHARACTER(len=256), PUBLIC :: ggmfile,ggmfile0 = "grid_geometry"        ! Head of time traces (gamma_x,Q_x)
  CHARACTER(len=256), PUBLIC :: prffile,prffile0 = "profiler"        ! Head of time traces (gamma_x,Q_x)
  CHARACTER(len=256), PUBLIC :: input_fname
  INTEGER, PUBLIC            :: fidres,fidmsp,fidfld,fidttr ! FID for output
  INTEGER, PUBLIC            :: fidmom,fidggm, fidprf

  PUBLIC :: diag_readinputs, diag_outputinputs, diagnose

CONTAINS

  SUBROUTINE diag_readinputs
    !    Read the input parameters

    USE basic, ONLY : lu_in, dt
    USE prec_const
    IMPLICIT NONE

    NAMELIST /DIAGNOSTICS/ dtsave_0d, dtsave_1d, dtsave_2d, dtsave_3d, dtsave_5d
    NAMELIST /DIAGNOSTICS/ write_doubleprecision, write_gamma, write_hf, write_phi
    NAMELIST /DIAGNOSTICS/ write_Na00, write_Napj, write_Sapj
    NAMELIST /DIAGNOSTICS/ write_dens, write_fvel, write_temp
    NAMELIST /DIAGNOSTICS/ diag_mode

    READ(lu_in,diagnostics)

    ! set nsave variables from dtsave ones (time unit to steps)
    nsave_0d = CEILING(dtsave_0d/dt)
    nsave_1d = CEILING(dtsave_1d/dt)
    nsave_2d = CEILING(dtsave_2d/dt)
    nsave_3d = CEILING(dtsave_3d/dt)
    nsave_5d = CEILING(dtsave_5d/dt)

  END SUBROUTINE diag_readinputs

  SUBROUTINE init_outfile(comm,file0,file)
    USE basic,           ONLY: speak, jobnum, basic_outputinputs
    USE grid,            ONLY: grid_outputinputs
    USE geometry,        ONLY: geometry_outputinputs
    USE model,           ONLY: model_outputinputs
    USE closure,         ONLY: closure_outputinputs
    USE species,         ONLY: species_outputinputs
    USE collision,       ONLY: coll_outputinputs
    USE initial,         ONLY: initial_outputinputs
    USE time_integration,ONLY: time_integration_outputinputs
    USE parallel,        ONLY: parallel_outputinputs
    USE units,           ONLY: units_outputinputs
    IMPLICIT NONE
    !input
    INTEGER,            INTENT(IN)    :: comm
    CHARACTER(len=256), INTENT(IN)    :: file0
    CHARACTER(len=256), INTENT(OUT)   :: file
    CHARACTER(len=256)                :: str
    INCLUDE 'srcinfo/srcinfo.h'
    ! Writing output filename
    WRITE(file,'(a,a1,i2.2,a3)') TRIM(file0)   ,'_',jobnum,'.h5'
    !                      1.1   Initial run
    CALL speak(TRIM(file)//' created',2)
    ! Add the code info and parameters to the file
    call h5write('codeinfo.h5', '/version', VERSION)
    call h5write('codeinfo.h5', '/branch', BRANCH)
    call h5write('codeinfo.h5', '/author', AUTHOR)
    call h5write('codeinfo.h5', '/execdate', EXECDATE)
    call h5write('codeinfo.h5', '/host', HOST)
    CALL            basic_outputinputs
    CALL             grid_outputinputs
    CALL         geometry_outputinputs
    CALL             diag_outputinputs
    CALL            model_outputinputs
    CALL          closure_outputinputs
    CALL          species_outputinputs
    CALL             coll_outputinputs
    CALL          initial_outputinputs
    CALL time_integration_outputinputs
    CALL         parallel_outputinputs
    CALL            units_outputinputs
    !  Save STDIN (input file) of this run
    IF(jobnum .LE. 99) THEN
       WRITE(str,'(a,i2.2)') "/files/STDIN.",jobnum
    ELSE
       WRITE(str,'(a,i3.2)') "/files/STDIN.",jobnum
    END IF
  END SUBROUTINE init_outfile

  SUBROUTINE diagnose(kstep)
    !   Diagnostics, writing simulation state to disk
    USE basic,           ONLY: lu_in, chrono_runt, cstep, time, display_h_min_s
    USE processing,      ONLY: pflux_x, hflux_x
    USE parallel,        ONLY: my_id
    USE units,           ONLY: WRITE_MW, pow_ref
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
    IF ((kstep .GE. 0) .AND. (MOD(cstep, nsave_0d) == 0) .AND. (my_id .EQ. 0)) THEN
      IF (WRITE_MW) THEN ! Display MW
        if(SIZE(pflux_x) .GT. 1) THEN
          WRITE(*,"(A,F8.2,A8,G10.3,A14,G10.3,A7)")&
          '|t = ', time,'| Qxi = ',hflux_x(1)*pow_ref,' [MW] | Qxe = ',hflux_x(2)*pow_Ref,' [MW] |'
        else
          WRITE(*,"(A,F8.2,A8,G10.3,A8,G10.3,A)")&
          '|t = ', time,'| Qxi = ',hflux_x(1)*pow_ref,' [MW]|'
        endif
      ELSE ! Display normalized code values
          if(SIZE(pflux_x) .GT. 1) THEN
            WRITE(*,"(A,F8.2,A8,G10.3,A8,G10.3,A8,G10.3,A8,G10.3,A)")&
            '|t = ', time,'| Pxi = ',pflux_x(1),'| Qxi = ',hflux_x(1),'| Pxe = ',pflux_x(2),'| Qxe = ',hflux_x(2),'|'
          else
            WRITE(*,"(A,F8.2,A8,G10.3,A8,G10.3,A)")&
            '|t = ', time,'| Pxi = ',pflux_x(1),'| Qxi = ',hflux_x(1),'|'
          endif
      ENDIF
    ENDIF

  END SUBROUTINE diagnose
  
  SUBROUTINE diagnose_full(kstep)
    USE basic,           ONLY: speak,chrono_runt,&
                               cstep,iframe0d,iframe2d,iframe3d,iframe5d,crashed
    USE grid,            ONLY: &
      parray_full,jarray_full, kparray, &
      kyarray_full,kxarray_full,zarray_full, ngz_o2, total_nz, local_nz, ieven,&
      local_Nky, total_nky, local_nkx, total_nkx
    USE geometry, ONLY: gxx, gxy, gxz, gyy, gyz, gzz, &
                        hatR, hatZ, hatB, dBdx, dBdy, dBdz, Jacobian, gradz_coeff
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
      CALL init_outfile(comm0,resfile0,resfile)
      ! Grid info
      CALL h5write("grid_geometry.h5", "/grid/local_nz", local_nz)
      CALL h5write("grid_geometry.h5", "/grid/total_nky", total_nky)
      CALL h5write("grid_geometry.h5", "/grid/local_Nky", local_Nky)
      CALL h5write("grid_geometry.h5", "/grid/total_nkx", total_nkx)
      CALL h5write("grid_geometry.h5", "/grid/coordkx", kxarray_full)
      CALL h5write("grid_geometry.h5", "/grid/coordky", kyarray_full)
      CALL h5write("grid_geometry.h5", "/grid/coordz", zarray_full)
      CALL h5write("grid_geometry.h5", "/grid/coordp", parray_full)
      CALL h5write("grid_geometry.h5", "/grid/coordj", jarray_full)
      CALL gather_xyz_real(kparray(1:local_Nky,1:local_Nkx,1:local_nz,ieven),kp_full,local_nky,total_nky,total_nkx,local_nz,total_nz)
      CALL h5write("grid_geometry.h5", "/grid/coordkp", kp_full)
      ! Metric info
      CALL gather_z(gxx((1+ngz_o2):(local_nz+ngz_o2),ieven),Az_full,local_nz,total_nz)
      CALL h5write("grid_geometry.h5", "/metric/gxx", Az_full)
      CALL gather_z(gxy((1+ngz_o2):(local_nz+ngz_o2),ieven),Az_full,local_nz,total_nz)
      CALL h5write("grid_geometry.h5", "/metric/gxy", Az_full)
      CALL gather_z(gxz((1+ngz_o2):(local_nz+ngz_o2),ieven),Az_full,local_nz,total_nz)
      CALL h5write("grid_geometry.h5", "/metric/gxz", Az_full)
      CALL gather_z(gyy((1+ngz_o2):(local_nz+ngz_o2),ieven),Az_full,local_nz,total_nz)
      CALL h5write("grid_geometry.h5", "/metric/gyy", Az_full)
      CALL gather_z(gyz((1+ngz_o2):(local_nz+ngz_o2),ieven),Az_full,local_nz,total_nz)
      CALL h5write("grid_geometry.h5", "/metric/gyz", Az_full)
      CALL gather_z(gzz((1+ngz_o2):(local_nz+ngz_o2),ieven),Az_full,local_nz,total_nz)
      CALL h5write("grid_geometry.h5", "/metric/gzz", Az_full)
      CALL gather_z(hatR((1+ngz_o2):(local_nz+ngz_o2),ieven),Az_full,local_nz,total_nz)
      CALL h5write("grid_geometry.h5", "/metric/hatR", Az_full)
      CALL gather_z(hatZ((1+ngz_o2):(local_nz+ngz_o2),ieven),Az_full,local_nz,total_nz)
      CALL h5write("grid_geometry.h5", "/metric/hatZ", Az_full)
      CALL gather_z(hatB((1+ngz_o2):(local_nz+ngz_o2),ieven),Az_full,local_nz,total_nz)
      CALL h5write("grid_geometry.h5", "/metric/hatB", Az_full)
      CALL gather_z(dBdx((1+ngz_o2):(local_nz+ngz_o2),ieven),Az_full,local_nz,total_nz)
      CALL h5write("grid_geometry.h5", "/metric/dBdx", Az_full)
      CALL gather_z(dBdy((1+ngz_o2):(local_nz+ngz_o2),ieven),Az_full,local_nz,total_nz)
      CALL h5write("grid_geometry.h5", "/metric/dBdy", Az_full)
      CALL gather_z(dBdz((1+ngz_o2):(local_nz+ngz_o2),ieven),Az_full,local_nz,total_nz)
      CALL h5write("grid_geometry.h5", "/metric/dBdz", Az_full)
      CALL gather_z(Jacobian((1+ngz_o2):(local_nz+ngz_o2),ieven),Az_full,local_nz,total_nz)
      CALL h5write("grid_geometry.h5", "/metric/Jacobian", Az_full)
      CALL gather_z(gradz_coeff((1+ngz_o2):(local_nz+ngz_o2),ieven),Az_full,local_nz,total_nz)
      CALL h5write("grid_geometry.h5", "/metric/gradz_coeff", Az_full)
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
      ! IF (nsave_2d .GT. 0) THEN
      !   IF (MOD(cstep, nsave_2d) == 0) THEN
      !     CALL diagnose_2d
      !   ENDIF
      ! ENDIF
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
  
  !!-------------- Auxiliary routines -----------------!!
  SUBROUTINE diagnose_0d
    USE basic
    USE prec_const
    USE processing, ONLY: pflux_x, hflux_x
    USE parallel, ONLY: my_id
    IMPLICIT NONE
    ! Write in a text file the time and the fluxes
    IF (my_id .EQ. 0) THEN
      WRITE(1,*) time, pflux_x(1), hflux_x(1)
    ENDIF
    ! write in a text file the profiler data
    IF (my_id .EQ. 0) THEN
      WRITE(2,*) time, chrono_mrhs%ttot, chrono_advf%ttot, chrono_clos%ttot, chrono_ghst%ttot, &
        chrono_coll%ttot, chrono_pois%ttot, chrono_sapj%ttot, chrono_chck%ttot, chrono_diag%ttot, &
        chrono_grad%ttot, chrono_napj%ttot, chrono_ExBs%ttot, chrono_step%ttot
    ENDIF
  END SUBROUTINE diagnose_0d
  
!   SUBROUTINE diagnose_2d
!     USE prec_const
!     USE basic
!     USE fields, ONLY: phi, psi, moments
!     USE array,  ONLY: Napjz,dens,upar,temp
!     USE grid, ONLY: CONTAINSp0, ip0,ij0, local_na, total_na,&
!          total_np, total_nj, total_nky, total_nkx, total_nz, &
!          local_np, local_nky, local_nkx, local_nz, &
!          ngz_o2, iz_obmp
!     USE time_integration, ONLY: updatetlevel
!     USE prec_const
!     USE processing, ONLY: compute_fluid_moments, compute_Napjz_spectrum
!     USE model,      ONLY: EM
!     USE parallel,   ONLY: manual_3D_bcast
! #ifdef TEST_SVD
!     USE CLA, ONLY: Sf
! #endif
!     IMPLICIT NONE
!     COMPLEX(xp), DIMENSION(local_na,local_nky,local_nkx,local_nz) :: Na00_
!     iframe2d=iframe2d+1
!     CALL attach(fidres,"/data/var2d/" , "frames", iframe2d)
!     CALL append(fidres,"/data/var2d/time", REAL(time,dp), ionode=0)
!     CALL append(fidres,"/data/var2d/cstep",REAL(cstep,dp),ionode=0)
! #ifdef TEST_SVD
!     WRITE(dset_name, "(A, '/', i6.6)") "/data/var2d/sv_ky_pj/", iframe2d
!     CALL putarr(fidres, dset_name, Sf, ionode=0)
! #endif
!     ! Write current EM fields
!     IF (write_phi)        CALL write_field2d_kykx_obmp(phi (:,:,(1+ngz_o2):(local_nz+ngz_o2)), 'phi_obmp')
!     IF (write_phi.AND.EM) CALL write_field2d_kykx_obmp(psi (:,:,(1+ngz_o2):(local_nz+ngz_o2)), 'psi_obmp')
!     IF (write_Na00) THEN
!       CALL compute_Napjz_spectrum
!       IF (CONTAINSp0) THEN
!         ! gyrocenter density
!         Na00_    = moments(:,ip0,ij0,:,:,(1+ngz_o2):(local_nz+ngz_o2),updatetlevel)
!       ELSE
!         Na00_    = 0._xp
!       ENDIF
!       CALL write_field2da_kykx_obmp(Na00_, 'Na00_obmp')
!       ! <<Napj>x>y spectrum
!       CALL write_field2da_pj_obmp(Napjz, 'Napj_obmp')
!     ENDIF
!     !! Fuid moments
!     IF (write_dens .OR. write_fvel .OR. write_temp) &
!     CALL compute_fluid_moments
!     IF (write_dens) CALL write_field2da_kykx_obmp(dens, 'dens_obmp')
!     IF (write_fvel) CALL write_field2da_kykx_obmp(upar, 'upar_obmp')
!     IF (write_temp) CALL write_field2da_kykx_obmp(temp, 'temp_obmp')
!     CONTAINS
!       SUBROUTINE write_field2d_kykx_obmp(field, text)
!         USE parallel, ONLY : gather_xyz, num_procs
!         IMPLICIT NONE
!         COMPLEX(xp), DIMENSION(:,:,:), INTENT(IN) :: field
!         CHARACTER(*), INTENT(IN) :: text
!         COMPLEX(xp), DIMENSION(total_nky,total_nkx,total_nz) :: field_full
!         CHARACTER(50) :: dset_name
!         field_full = 0;
!         WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var2d", TRIM(text), iframe2d
!         IF (num_procs .EQ. 1) THEN ! no data distribution
!           CALL putarr(fidres, dset_name, field(:,:,iz_obmp), ionode=0)
!         ELSE
!           CALL gather_xyz(field,field_full,local_nky,total_nky,total_nkx,local_nz,total_nz)
!           CALL putarr(fidres, dset_name, field_full(:,:,iz_obmp), ionode=0)
!         ENDIF
!         END SUBROUTINE write_field2d_kykx_obmp
  
!         SUBROUTINE write_field2da_kykx_obmp(field, text)
!           USE parallel, ONLY : gather_xyz, num_procs
!           IMPLICIT NONE
!           COMPLEX(xp), DIMENSION(:,:,:,:), INTENT(IN) :: field
!           CHARACTER(*), INTENT(IN) :: text
!           COMPLEX(xp), DIMENSION(total_na,total_nky,total_nkx,total_nz) :: field_full
!           COMPLEX(xp), DIMENSION(local_nky,total_nkx,local_nz) :: buff_local
!           COMPLEX(xp), DIMENSION(total_nky,total_nkx,total_nz) :: buff_full
!           CHARACTER(50) :: dset_name
!           INTEGER :: ia
!           field_full = 0;
!           WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var2d", TRIM(text), iframe2d
!           IF (num_procs .EQ. 1) THEN ! no data distribution
!             CALL putarr(fidres, dset_name, field(:,:,:,iz_obmp), ionode=0)
!           ELSE
!             DO ia = 1,total_Na
!               buff_local = field(ia,:,:,:)
!               CALL gather_xyz(buff_local,buff_full,local_nky,total_nky,total_nkx,local_nz,total_nz)
!               field_full(ia,:,:,:) = buff_full
!             ENDDO
!             CALL putarr(fidres, dset_name, field_full(:,:,:,iz_obmp), ionode=0)
!           ENDIF
!           END SUBROUTINE write_field2da_kykx_obmp
  
!       SUBROUTINE write_field2da_pj_obmp(field, text)
!         USE parallel, ONLY : gather_pjz, num_procs
!         IMPLICIT NONE
!         COMPLEX(xp), DIMENSION(:,:,:,:), INTENT(IN) :: field
!         CHARACTER(*), INTENT(IN) :: text
!         COMPLEX(xp), DIMENSION(total_na,total_np,total_nj,total_nz) :: field_full
!         CHARACTER(LEN=50) :: dset_name
!         field_full = 0;
!         WRITE(dset_name, "(A, '/', A, '/', i6.6)") "/data/var2d", TRIM(text), iframe2d
!         IF (num_procs .EQ. 1) THEN ! no data distribution
!           CALL putarr(fidres, dset_name, field(:,:,:,iz_obmp), ionode=0)
!         ELSE
!           CALL gather_pjz(field,field_full,total_na,local_np,total_np,total_nj,local_nz,total_nz)
!           CALL putarr(fidres, dset_name, field_full(:,:,:,iz_obmp), ionode=0)
!         ENDIF
!         CALL attach(fidres, dset_name, "time", time)
!       END SUBROUTINE write_field2da_pj_obmp
!   END SUBROUTINE diagnose_2d
  
  SUBROUTINE diagnose_3d
    USE basic
    USE fields, ONLY: phi, psi, moments
    USE array,  ONLY: Napjz,dens,upar,uper,Tpar,Tper,temp
    USE grid, ONLY: CONTAINSp0, ip0,ij0, local_na, total_na,&
                    total_np, total_nj, total_nky, total_nkx, total_nz, &
                    local_np, local_nky, local_nkx, local_nz, &
                    ngz_o2
    USE time_integration, ONLY: updatetlevel
    USE prec_const
    USE processing, ONLY: compute_fluid_moments, compute_Napjz_spectrum
    USE model,      ONLY: EM
    USE parallel,   ONLY: manual_3D_bcast
    IMPLICIT NONE
    COMPLEX(xp), DIMENSION(local_na,local_nky,local_nkx,local_nz) :: Na00_
    CHARACTER(50) :: dset_name, filename
    iframe3d=iframe3d+1
    WRITE(filename, "(A, i6.6, A)") "diag3d_", iframe3d, ".h5"    
    ! add current time, cstep and frame
    CALL h5write(filename, "/time", time)
    CALL h5write(filename, "/cstep", cstep)
    CALL h5write(filename, "/frame", iframe3d)
    ! Write current EM fields
    IF (write_phi)        CALL write_field3d_kykxz(phi (:,:,(1+ngz_o2):(local_nz+ngz_o2)), 'phi')
    IF (write_phi.AND.EM) CALL write_field3d_kykxz(psi (:,:,(1+ngz_o2):(local_nz+ngz_o2)), 'psi')
    IF (write_Na00) THEN
      CALL compute_Napjz_spectrum
      IF (CONTAINSp0) THEN
        ! gyrocenter density
        Na00_    = moments(:,ip0,ij0,:,:,(1+ngz_o2):(local_nz+ngz_o2),updatetlevel)
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
        CHARACTER(50) :: dset_name, filename
        WRITE(filename, "(A, i6.6, A)") "diag3d_", iframe3d, ".h5"
        WRITE(dset_name, "(A, A)") "/", TRIM(text)
        field_full = 0;
        IF (num_procs .EQ. 1) THEN ! no data distribution
          CALL h5write(filename, dset_name, field)
        ELSE
          CALL gather_xyz(field,field_full,local_nky,total_nky,total_nkx,local_nz,total_nz)
          CALL h5write(filename, dset_name, field_full)
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
          CHARACTER(50) :: dset_name, filename
          INTEGER :: ia
          WRITE(filename, "(A, i6.6, A)") "diag3d_", iframe3d, ".h5"
          WRITE(dset_name, "(A, A)") "/", TRIM(text)
          field_full = 0;
          IF (num_procs .EQ. 1) THEN ! no data distribution
            CALL h5write(filename, dset_name, field)
          ELSE
            DO ia = 1,total_Na
              buff_local = field(ia,:,:,:)
              CALL gather_xyz(buff_local,buff_full,local_nky,total_nky,total_nkx,local_nz,total_nz)
              field_full(ia,:,:,:) = buff_full
            ENDDO
            CALL h5write(filename, dset_name, field_full)
          ENDIF
          END SUBROUTINE write_field3da_kykxz
  
      SUBROUTINE write_field3da_pjz(field, text)
        USE parallel, ONLY : gather_pjz, num_procs
        IMPLICIT NONE
        COMPLEX(xp), DIMENSION(:,:,:,:), INTENT(IN) :: field
        CHARACTER(*), INTENT(IN) :: text
        COMPLEX(xp), DIMENSION(total_na,total_np,total_nj,total_nz) :: field_full
        CHARACTER(LEN=50) :: dset_name
        CHARACTER(LEN=50) :: filename
        field_full = 0;
        WRITE(filename, "(A, i6.6, A)") "diag3d_", iframe3d, ".h5"
        WRITE(dset_name, "(A, A)") "/", TRIM(text)
        IF (num_procs .EQ. 1) THEN ! no data distribution
          CALL h5write(filename, dset_name, field)
        ELSE
          CALL gather_pjz(field,field_full,total_na,local_np,total_np,total_nj,local_nz,total_nz)
          CALL h5write(filename, dset_name, field_full)
        ENDIF
      END SUBROUTINE write_field3da_pjz
  
  END SUBROUTINE diagnose_3d
  
  SUBROUTINE diagnose_5d
  
    USE basic,  ONLY: time, iframe5d,cstep
    USE fields, ONLY: moments
    USE grid,   ONLY:total_np, total_nj, total_nky, total_nkx, total_nz, &
                     local_np, local_nj, local_nky, local_nkx, local_nz, &
                     ngp, ngj, ngz, ngz_o2, total_na
    USE prec_const, ONLY: xp, dp
    IMPLICIT NONE
    iframe5d=iframe5d+1
    IF (write_Napj) THEN
    CALL write_field5d(moments)
    ENDIF
  
    CONTAINS
  
    SUBROUTINE write_field5d(field)
      USE basic,            ONLY: jobnum, dt
      USE parallel,         ONLY: gather_pjxyz, num_procs
      USE prec_const,       ONLY: xp
      USE time_integration, ONLY: ntimelevel
      IMPLICIT NONE
      COMPLEX(xp), DIMENSION(total_na,local_np+ngp,local_nj+ngj,local_Nky,local_Nkx,local_nz+ngz,ntimelevel), INTENT(IN) :: field
      COMPLEX(xp), DIMENSION(total_na,local_np,local_nj,local_nky,local_nkx,local_nz) :: field_sub
      COMPLEX(xp), DIMENSION(total_na,total_np,total_nj,total_nky,total_nkx,total_nz) :: field_full
      CHARACTER(LEN=50) :: filename
      field_sub  = field(1:total_na,(1+ngp/2):(local_np+ngp/2),(1+ngj/2):(local_nj+ngj/2),&
                            1:local_nky,1:local_nkx,(1+ngz_o2):(local_nz+ngz_o2),1)
      field_full = 0;
      WRITE(filename, "(A, i6.6, A)") "moments_", iframe5d, ".h5"
      IF (num_procs .GT. 1) THEN
        CALL gather_pjxyz(field_sub,field_full,total_na,local_np,total_np,total_nj,local_nky,total_nky,total_nkx,local_nz,total_nz)
        CALL h5write(filename, "/moments", field_full)
      ELSE
        CALL h5write(filename, "/moments", field_sub)
      ENDIF
      CALL h5write(filename, "/time", time)
      CALL h5write(filename, "/jobnum", jobnum)
      CALL h5write(filename, "/dt", dt)
      CALL h5write(filename, "/cstep", cstep)
      CALL h5write(filename, "/iframe5d", iframe5d)
    END SUBROUTINE write_field5d
  END SUBROUTINE diagnose_5d  

  SUBROUTINE diag_outputinputs
    !
    !    Write the input parameters to the results_xx.h5 file
    !
    USE h5fortran
    IMPLICIT NONE
    CALL h5write("outputinput.h5", "/diagnostics/write_doubleprecision", write_doubleprecision)
    CALL h5write("outputinput.h5", "/diagnostics/nsave_0d", nsave_0d)
    CALL h5write("outputinput.h5", "/diagnostics/nsave_1d", nsave_1d)
    CALL h5write("outputinput.h5", "/diagnostics/nsave_2d", nsave_2d)
    CALL h5write("outputinput.h5", "/diagnostics/nsave_3d", nsave_3d)
    CALL h5write("outputinput.h5", "/diagnostics/nsave_5d", nsave_5d)
    CALL h5write("outputinput.h5", "/diagnostics/write_gamma", write_gamma)
    CALL h5write("outputinput.h5", "/diagnostics/write_hf", write_hf)
    CALL h5write("outputinput.h5", "/diagnostics/write_phi", write_phi)
    CALL h5write("outputinput.h5", "/diagnostics/write_Na00", write_Na00)
    CALL h5write("outputinput.h5", "/diagnostics/write_Napj", write_Napj)
    CALL h5write("outputinput.h5", "/diagnostics/write_Sapj", write_Sapj)
    CALL h5write("outputinput.h5", "/diagnostics/write_dens", write_dens)
    CALL h5write("outputinput.h5", "/diagnostics/write_fvel", write_fvel)
    CALL h5write("outputinput.h5", "/diagnostics/write_temp", write_temp)
  END SUBROUTINE diag_outputinputs

END MODULE diagnostics
