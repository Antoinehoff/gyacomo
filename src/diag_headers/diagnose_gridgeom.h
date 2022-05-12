SUBROUTINE diagnose_gridgeom(kstep)
  USE basic
  USE grid
  USE diagnostics_par
  USE futils, ONLY: putarr, creatg, putarrnd, closef
  USE time_integration
  USE prec_const
  USE geometry
  USE array
  IMPLICIT NONE
  INTEGER, INTENT(in) :: kstep
  INTEGER, parameter  :: BUFSIZE = 2
  INTEGER :: rank = 0
  INTEGER :: dims(1) = (/0/)
  IF (kstep .EQ. 0) THEN
    ! Grid info
    CALL init_outfile(comm0,ggmfile0,ggmfile,fidggm)
    CALL creatg(fidggm, "/data/grid", "Grid data")
    CALL putarr(fidggm, "/data/grid/coordkx",   kxarray_full,  "kx*rho_s0", ionode=0)
    CALL putarr(fidggm, "/data/grid/coordky",   kyarray_full,  "ky*rho_s0", ionode=0)
    CALL putarr(fidggm, "/data/grid/coordz",    zarray_full,   "z/R", ionode=0)
    CALL putarr(fidggm, "/data/grid/coordp_e" , parray_e_full, "p_e", ionode=0)
    CALL putarr(fidggm, "/data/grid/coordj_e" , jarray_e_full, "j_e", ionode=0)
    CALL putarr(fidggm, "/data/grid/coordp_i" , parray_i_full, "p_i", ionode=0)
    CALL putarr(fidggm, "/data/grid/coordj_i" , jarray_i_full, "j_i", ionode=0)

    ! Metric info
    CALL   creatg(fidggm, "/data/metric", "Metric data")
    CALL putarrnd(fidggm, "/data/metric/gxx",            gxx(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidggm, "/data/metric/gxy",            gxy(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidggm, "/data/metric/gyy",            gyy(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidggm, "/data/metric/gyz",            gyz(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidggm, "/data/metric/gzz",            gzz(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidggm, "/data/metric/hatR",          hatR(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidggm, "/data/metric/hatZ",          hatZ(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidggm, "/data/metric/hatB",          hatB(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidggm, "/data/metric/gradxB",      gradxB(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidggm, "/data/metric/gradyB",      gradyB(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidggm, "/data/metric/gradzB",      gradzB(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidggm, "/data/metric/Jacobian",    Jacobian(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidggm, "/data/metric/gradz_coeff", gradz_coeff(izs:ize,0:1), (/1, 1, 1/))
    CALL putarrnd(fidggm, "/data/metric/Ckxky",       Ckxky(ikys:ikye,ikxs:ikxe,izs:ize,0:1), (/1, 1, 3/))
    CALL putarrnd(fidggm, "/data/metric/kernel_i",    kernel_i(ijs_i:ije_i,ikys:ikye,ikxs:ikxe,izs:ize,0:1), (/ 1, 2, 4/))
  ENDIF
  IF (kstep .EQ. -1) THEN
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    CALL closef(fidggm)
  ENDIF

END SUBROUTINE diagnose_gridgeom
