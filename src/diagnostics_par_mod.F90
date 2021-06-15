MODULE diagnostics_par
  !   Module for diagnostic parameters

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PUBLIC, PROTECTED :: write_doubleprecision = .FALSE.
  LOGICAL, PUBLIC, PROTECTED :: write_gamma
  LOGICAL, PUBLIC, PROTECTED :: write_phi,  write_Na00
  LOGICAL, PUBLIC, PROTECTED :: write_Napj, write_Sapj
  LOGICAL, PUBLIC, PROTECTED :: write_dens, write_temp

  INTEGER, PUBLIC, PROTECTED :: nsave_0d , nsave_1d , nsave_2d , nsave_5d, nsave_cp

  !  HDF5 file
  CHARACTER(len=256), PUBLIC :: resfile0 = "results"   ! Head of main result file name
  CHARACTER(len=256), PUBLIC :: resfile                ! Main result file
  INTEGER, PUBLIC            :: job2load               ! jobnum of the checkpoint to load
  INTEGER, PUBLIC            :: fidres                 ! FID for resfile
  CHARACTER(len=256), PUBLIC :: rstfile0 = "restart"   ! Head of restart file name
  CHARACTER(len=256), PUBLIC :: rstfile                ! Full restart file
  INTEGER, PUBLIC            :: fidrst                 ! FID for restart file

  PUBLIC :: output_par_readinputs, output_par_outputinputs

CONTAINS


  SUBROUTINE output_par_readinputs
    !    Read the input parameters

    USE basic, ONLY : lu_in
    USE prec_const
    IMPLICIT NONE

    NAMELIST /OUTPUT_PAR/ nsave_0d , nsave_1d , nsave_2d , nsave_5d, nsave_cp
    NAMELIST /OUTPUT_PAR/ write_doubleprecision, write_gamma, write_phi
    NAMELIST /OUTPUT_PAR/ write_Na00, write_Napj, write_Sapj
    NAMELIST /OUTPUT_PAR/ write_dens, write_temp
    NAMELIST /OUTPUT_PAR/ resfile0, rstfile0, job2load

    READ(lu_in,output_par)

  END SUBROUTINE output_par_readinputs


  SUBROUTINE output_par_outputinputs(fidres, str)
    !
    !    Write the input parameters to the results_xx.h5 file
    !
    USE prec_const
    USE futils, ONLY: attach
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fidres
    CHARACTER(len=256), INTENT(in) :: str

    CALL attach(fidres, TRIM(str), "write_doubleprecision", write_doubleprecision)
    CALL attach(fidres, TRIM(str), "nsave_0d", nsave_0d)
    CALL attach(fidres, TRIM(str), "nsave_1d", nsave_1d)
    CALL attach(fidres, TRIM(str), "nsave_2d", nsave_2d)
    CALL attach(fidres, TRIM(str), "nsave_5d", nsave_5d)
    CALL attach(fidres, TRIM(str), "resfile0", resfile0)

  END SUBROUTINE output_par_outputinputs


END MODULE diagnostics_par
