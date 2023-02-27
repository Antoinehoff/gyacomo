subroutine auxval
  !   Set auxiliary values, at beginning of simulation

  USE basic
  USE grid
  USE array
  USE model
  USE fourier, ONLY: init_grid_distr_and_plans
  use prec_const
  USE numerics
  USE geometry
  USE parallel, ONLY: init_parallel_var
  IMPLICIT NONE

  INTEGER :: i_
  IF (my_id .EQ. 0) WRITE(*,*) '=== Set auxiliary values ==='

  IF (LINEARITY .NE. 'linear') THEN
    IF (my_id .EQ. 0) write(*,*) 'FFTW3 y-grid distribution'
    CALL init_grid_distr_and_plans(Nx,Ny)
  ELSE
    CALL init_1Dgrid_distr
    IF (my_id .EQ. 0) write(*,*) 'Manual y-grid distribution'
  ENDIF
  ! Init the grids
  CALL set_pgrid ! parallel kin (MPI distributed)

  CALL set_jgrid ! perp kin

  CALL set_kxgrid(shear) ! radial modes (MPI distributed by FFTW)

  CALL set_kygrid ! azymuthal modes

  CALL set_zgrid  ! field aligned angle
  IF ((my_id .EQ. 0) .AND. SG) WRITE(*,*) '--2 staggered z grids--'

  CALL memory ! Allocate memory for global arrays

  CALL init_parallel_var

  CALL eval_magnetic_geometry ! precompute coeff for lin equation

  CALL compute_lin_coeff ! precompute coeff for lin equation and geometry

  CALL evaluate_kernels ! precompute the kernels

  CALL evaluate_EM_op ! compute inverse of poisson and ampere operators

  IF ( LINEARITY .NE. 'linear' ) THEN;
    CALL build_dnjs_table ! precompute the Laguerre nonlin product coeffs
  ENDIF

  CALL build_dv4Hp_table ! precompute the hermite fourth derivative table

  !! Display parallel settings
  DO i_ = 0,num_procs-1
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    IF (my_id .EQ. i_) THEN
      IF (my_id .EQ. 0) WRITE(*,*) ''
      IF (my_id .EQ. 0) WRITE(*,*) '--------- Parallel environement ----------'
      IF (my_id .EQ. 0) WRITE(*,'(A12,I3)') 'n_procs ', num_procs
      IF (my_id .EQ. 0) WRITE(*,'(A12,I3,A14,I3,A14,I3)') 'num_procs_p   = ', num_procs_p, ', num_procs_ky   = ', num_procs_ky, ', num_procs_z   = ', num_procs_z
      IF (my_id .EQ. 0) WRITE(*,*) ''
      WRITE(*,'(A9,I3,A10,I3,A10,I3,A9,I3)')&
       'my_id  = ', my_id, ', rank_p  = ', rank_p, ', rank_ky  = ', rank_ky,', rank_z  = ', rank_z
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
      IF (my_id .NE. num_procs-1) WRITE (*,*) ''
      IF (my_id .EQ. num_procs-1) WRITE(*,*) '------------------------------------------'
    ENDIF
  ENDDO
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)

  IF((my_id.EQ.0) .AND. (CLOS .EQ. 1)) THEN
  IF(KIN_E) &
  write(*,*) 'Closure = 1 -> Maximal Nepj degree is min(Pmaxe,2*Jmaxe+1): De = ', dmaxi
  write(*,*) 'Closure = 1 -> Maximal Nipj degree is min(Pmaxi,2*Jmaxi+1): Di = ', dmaxi
  ENDIF

END SUBROUTINE auxval
