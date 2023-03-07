MODULE diagnostics_par
  !   Module for diagnostic parameters

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PUBLIC, PROTECTED :: write_doubleprecision = .FALSE.
  LOGICAL, PUBLIC, PROTECTED :: write_gamma, write_hf  ! output particle transport and heat flux
  LOGICAL, PUBLIC, PROTECTED :: write_phi,  write_Na00
  LOGICAL, PUBLIC, PROTECTED :: write_Napj, write_Sapj
  LOGICAL, PUBLIC, PROTECTED :: write_dens, write_fvel, write_temp

  INTEGER, PUBLIC, PROTECTED :: nsave_0d, nsave_1d, nsave_2d, nsave_3d, nsave_5d ! save data every n step
  REAL,    PUBLIC, PROTECTED :: dtsave_0d, dtsave_1d, dtsave_2d, dtsave_3d, dtsave_5d ! save data every dt time unit

  !  HDF5 file
  CHARACTER(len=256), PUBLIC :: resfile,resfile0 = "outputs"            ! Head of main result file name
  CHARACTER(len=256), PUBLIC :: momfile,momfile0 = "moments"   ! Head of the moment spectrum file (N_a(p,j,z))
  CHARACTER(len=256), PUBLIC :: mspfile,mspfile0 = "moments_spectrum"   ! Head of the moment spectrum file (N_a(p,j,z))
  CHARACTER(len=256), PUBLIC :: fldfile,fldfile0 = "fields"             ! Head of field (phi,A)
  CHARACTER(len=256), PUBLIC :: ttrfile,ttrfile0 = "time_traces"        ! Head of time traces (gamma_x,Q_x)
  CHARACTER(len=256), PUBLIC :: ggmfile,ggmfile0 = "grid_geometry"        ! Head of time traces (gamma_x,Q_x)
  CHARACTER(len=256), PUBLIC :: prffile,prffile0 = "profiler"        ! Head of time traces (gamma_x,Q_x)
  CHARACTER(len=256), PUBLIC :: input_fname
  CHARACTER(len=256), PUBLIC :: rstfile                     ! restart result file
  INTEGER, PUBLIC            :: fidres,fidmsp,fidfld,fidttr ! FID for output
  INTEGER, PUBLIC            :: fidmom,fidggm, fidprf
  INTEGER, PUBLIC            :: fidrst                      ! FID for restart file

  PUBLIC :: diag_par_readinputs, diag_par_outputinputs

CONTAINS


  SUBROUTINE diag_par_readinputs
    !    Read the input parameters

    USE basic, ONLY : lu_in, dt
    USE prec_const
    IMPLICIT NONE

    NAMELIST /OUTPUT_PAR/ dtsave_0d, dtsave_1d, dtsave_2d, dtsave_3d, dtsave_5d
    NAMELIST /OUTPUT_PAR/ write_doubleprecision, write_gamma, write_hf, write_phi
    NAMELIST /OUTPUT_PAR/ write_Na00, write_Napj, write_Sapj
    NAMELIST /OUTPUT_PAR/ write_dens, write_fvel, write_temp

    READ(lu_in,output_par)

    ! set nsave variables from dtsave ones (time unit to steps)
    nsave_0d = CEILING(dtsave_0d/dt)
    nsave_1d = CEILING(dtsave_1d/dt)
    nsave_2d = CEILING(dtsave_2d/dt)
    nsave_3d = CEILING(dtsave_3d/dt)
    nsave_5d = CEILING(dtsave_5d/dt)

  END SUBROUTINE diag_par_readinputs


  SUBROUTINE diag_par_outputinputs(fid)
    !
    !    Write the input parameters to the results_xx.h5 file
    !
    USE prec_const
    USE futils, ONLY: attach, creatd
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fid
    CHARACTER(len=256)  :: str
    WRITE(str,'(a)') '/data/input/diag_par'
    CALL creatd(fidres, 0,(/0/),TRIM(str),'Diagnostics Parameters Input')
    CALL attach(fid, TRIM(str), "write_doubleprecision", write_doubleprecision)
    CALL attach(fid, TRIM(str), "nsave_0d", nsave_0d)
    CALL attach(fid, TRIM(str), "nsave_1d", nsave_1d)
    CALL attach(fid, TRIM(str), "nsave_2d", nsave_2d)
    CALL attach(fid, TRIM(str), "nsave_3d", nsave_3d)
    CALL attach(fid, TRIM(str), "nsave_5d", nsave_5d)
    CALL attach(fid, TRIM(str), "write_gamma", write_gamma)
    CALL attach(fid, TRIM(str),    "write_hf",    write_hf)
    CALL attach(fid, TRIM(str),   "write_phi",   write_phi)
    CALL attach(fid, TRIM(str),  "write_Na00",  write_Na00)
    CALL attach(fid, TRIM(str),  "write_Napj",  write_Napj)
    CALL attach(fid, TRIM(str),  "write_Sapj",  write_Sapj)
    CALL attach(fid, TRIM(str),  "write_dens",  write_dens)
    CALL attach(fid, TRIM(str),  "write_fvel",  write_fvel)
    CALL attach(fid, TRIM(str),  "write_temp",  write_temp)

  END SUBROUTINE diag_par_outputinputs


END MODULE diagnostics_par
