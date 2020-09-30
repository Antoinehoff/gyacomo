MODULE diagnostics_par
  !   Module for diagnostic parameters

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PUBLIC, PROTECTED :: write_Ni00=.TRUE.
  LOGICAL, PUBLIC, PROTECTED :: write_moments=.TRUE.
  LOGICAL, PUBLIC, PROTECTED :: write_phi=.TRUE.
  LOGICAL, PUBLIC, PROTECTED :: write_non_lin=.TRUE.
  LOGICAL, PUBLIC, PROTECTED :: write_doubleprecision=.FALSE.

  INTEGER, PUBLIC, PROTECTED :: nsave_0d , nsave_1d , nsave_2d , nsave_5d
  INTEGER, PUBLIC, PROTECTED :: nsave_cp = 1e3

  !  HDF5 file
  CHARACTER(len=64), PUBLIC :: resfile0 = "results"   ! Head of main result file name
  CHARACTER(len=64), PUBLIC :: resfile                ! Main result file
  INTEGER, PUBLIC           :: fidres                 ! FID for resfile
  CHARACTER(len=64), PUBLIC :: rstfile0 = "restart"   ! Head of restart file name
  CHARACTER(len=64), PUBLIC :: rstfile                ! Full restart file
  INTEGER, PUBLIC           :: fidrst                 ! FID for restart file

  PUBLIC :: output_par_readinputs, output_par_outputinputs

CONTAINS


  SUBROUTINE output_par_readinputs
    !    Read the input parameters

    USE basic, ONLY : lu_in
    USE prec_const
    IMPLICIT NONE

    NAMELIST /OUTPUT_PAR/ nsave_0d , nsave_1d , nsave_2d , nsave_5d
    NAMELIST /OUTPUT_PAR/ write_Ni00, write_moments, write_phi, write_non_lin, write_doubleprecision
    NAMELIST /OUTPUT_PAR/ resfile0, rstfile0

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

    CALL attach(fidres, TRIM(str), "write_Ni00", write_Ni00)
    CALL attach(fidres, TRIM(str), "write_moments", write_moments)
    CALL attach(fidres, TRIM(str), "write_phi", write_phi)
    CALL attach(fidres, TRIM(str), "write_non_lin", write_non_lin)
    CALL attach(fidres, TRIM(str), "write_doubleprecision", write_doubleprecision)
    CALL attach(fidres, TRIM(str), "nsave_0d", nsave_0d)
    CALL attach(fidres, TRIM(str), "nsave_1d", nsave_1d)
    CALL attach(fidres, TRIM(str), "nsave_2d", nsave_2d)
    CALL attach(fidres, TRIM(str), "nsave_5d", nsave_5d)
    CALL attach(fidres, TRIM(str), "resfile0", resfile0)

  END SUBROUTINE output_par_outputinputs


END MODULE diagnostics_par
