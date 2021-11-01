subroutine auxval
  !   Set auxiliary values, at beginning of simulation

  USE basic
  USE grid
  USE array
  USE model
  USE fourier, ONLY: init_grid_distr_and_plans, alloc_local_1, alloc_local_2
  use prec_const
  USE numerics
  USE geometry
  IMPLICIT NONE

  INTEGER :: irows,irowe, irow, icol, i_
  IF (my_id .EQ. 0) WRITE(*,*) '=== Set auxiliary values ==='

  IF (NON_LIN) THEN
    CALL init_grid_distr_and_plans(Nx,Ny)
  ELSE
    CALL init_1Dgrid_distr
  ENDIF
  ! Init the grids
  CALL set_pgrid ! parallel kin (MPI distributed)

  CALL set_jgrid ! perp kin

  CALL set_kxgrid ! radial modes (MPI distributed by FFTW)

  CALL set_kygrid ! azymuthal modes

  CALL set_zgrid  ! field aligned angle

  CALL memory ! Allocate memory for global arrays

  CALL eval_magnetic_geometry ! precompute coeff for lin equation

  CALL compute_lin_coeff ! precompute coeff for lin equation and geometry

  CALL evaluate_kernels ! precompute the kernels

  CALL evaluate_poisson_op ! precompute the kernels

  IF ( NON_LIN ) THEN;
    CALL build_dnjs_table ! precompute the Laguerre nonlin product coeffs
  ENDIF

  !! Display parallel settings
  DO i_ = 0,num_procs-1
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    IF (my_id .EQ. i_) THEN
      IF (my_id .EQ. 0) WRITE(*,*) ''
      IF (my_id .EQ. 0) WRITE(*,*) '--------- Parallel environement ----------'
      IF (my_id .EQ. 0) WRITE(*,'(A9,I3,A10,I3,A10,I3)') 'n_procs= ', num_procs, ', num_procs_p   = ', num_procs_p, ', num_procs_kx   = ', num_procs_kx
      IF (my_id .EQ. 0) WRITE(*,*) ''
      WRITE(*,'(A9,I3,A10,I3,A10,I3)')&
       'my_id  = ', my_id, ', rank_p  = ', rank_p, ', rank_kx  = ', rank_kx
       WRITE(*,'(A22,I3,A11,I3)')&
       '              ips_e = ', ips_e, ', ipe_e  = ', ipe_e
       WRITE(*,'(A22,I3,A11,I3)')&
       '              ijs_e = ', ijs_e, ', ije_e  = ', ije_e
       WRITE(*,'(A22,I3,A11,I3)')&
       '              ips_i = ', ips_i, ', ipe_i  = ', ipe_i
       WRITE(*,'(A22,I3,A11,I3)')&
       '              ijs_i = ', ijs_i, ', ije_i  = ', ije_i
       WRITE(*,'(A22,I3,A11,I3)')&
       '              ikxs  = ', ikxs , ', ikxe   = ', ikxe
       WRITE(*,'(A22,I3,A11,I3)')&
       '              ikys  = ', ikys , ', ikye   = ', ikye
       WRITE(*,'(A22,I3,A11,I3)')&
       '              izs   = ', izs  , ', ize    = ', ize
       ! write(*,*) 'local kx =', kxarray
       ! write(*,*) 'local ky =', kyarray
       ! write(*,*) 'local iz =', izarray
      IF (my_id .NE. num_procs-1) WRITE (*,*) ''
      IF (my_id .EQ. num_procs-1) WRITE(*,*) '------------------------------------------'
    ENDIF
  ENDDO
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)

END SUBROUTINE auxval
