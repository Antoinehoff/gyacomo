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
  USE parallel, ONLY: init_parallel_var, my_id, num_procs, num_procs_p, num_procs_z, num_procs_ky, rank_p, rank_ky, rank_z
  USE processing, ONLY: init_process
  IMPLICIT NONE

  INTEGER :: i_, ierr
  IF (my_id .EQ. 0) WRITE(*,*) '=== Set auxiliary values ==='

  ! Init the grids
  CALL set_grids(shear,Npol) ! radial modes (MPI distributed by FFTW)

  CALL memory ! Allocate memory for global arrays

  CALL init_parallel_var(local_np,total_np,local_nky,total_nky,local_nz)

  CALL init_process

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
       '              ips = ', ips, ', ipe  = ', ipe
       WRITE(*,'(A22,I3,A11,I3)')&
       '              ijs = ', ijs, ', ije  = ', ije
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

  IF((CLOS .EQ. 1)) &
    CALL speak('Closure = 1 -> Maximal Napj degree is min(Pmax,2*Jmax+1): D = '// str(dmax))

END SUBROUTINE auxval
