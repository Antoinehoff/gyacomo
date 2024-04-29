subroutine auxval
  !   Set auxiliary values, at beginning of simulation
  USE, intrinsic :: iso_fortran_env, ONLY: OUTPUT_UNIT
  USE basic,          ONLY: str, speak, VERBOSE_LVL
  USE grid,           ONLY: local_np, local_np_offset, total_np, local_nj, local_nj_offset,&
                            local_nky, local_nky_offset, total_nky, local_nkx, local_nkx_offset,&
                            local_nz, local_nz_offset, init_grids_data, set_grids
  !USE array
  USE model,          ONLY: Na, EM, LINEARITY, N_HD, N_HDz, ExBrate
  USE fourier,        ONLY: init_grid_distr_and_plans
  use MPI,            ONLY: MPI_COMM_WORLD
  USE numerics,       ONLY: build_dnjs_table, build_dv4Hp_table, compute_lin_coeff, &
                            evaluate_EM_op, evaluate_kernels
  USE geometry,       ONLY: Npol, shear, eval_magnetic_geometry
  USE closure,        ONLY: set_closure_model, hierarchy_closure, dmax
  USE parallel,       ONLY: init_parallel_var, my_id, num_procs, &
    num_procs_p, num_procs_z, num_procs_ky, rank_p, rank_ky, rank_z
  USE processing,     ONLY: init_process
  USE ExB_shear_flow, ONLY: Setup_ExB_shear_flow
#ifdef TEST_SVD
  USE CLA, ONLY: init_CLA
#endif
  IMPLICIT NONE
  INTEGER :: i_, ierr
  ! Init the grids
  CALL init_grids_data(Na,EM,LINEARITY) 
  CALL set_grids(shear,ExBrate,Npol,LINEARITY,N_HD,N_HDz,EM,Na) 
  ! Allocate memory for global arrays
  CALL memory
  ! Initialize displacement and receive arrays
  CALL init_parallel_var(local_np,total_np,local_nky,total_nky,local_nz)
  ! Initialize heatflux variables
  CALL init_process
  ! precompute coeff for lin equation
  CALL eval_magnetic_geometry 
  ! precompute coeff for lin equation and geometry
  CALL compute_lin_coeff 
  ! precompute the kernels
  CALL evaluate_kernels 
  ! compute inverse of poisson and ampere operators
  CALL evaluate_EM_op
  ! precompute the Laguerre nonlin product coeffs
  IF ( LINEARITY .NE. 'linear' ) &
    CALL build_dnjs_table 
  ! precompute the hermite fourth derivative table
  CALL build_dv4Hp_table 
  ! set the closure scheme in use
  CALL set_closure_model   
  ! Setup ExB shear variables
  CALL Setup_ExB_shear_flow(ExBrate)
#ifdef TEST_SVD
  ! If we want to test SVD decomposition etc.
  CALL init_CLA(local_nky,local_np*local_nj)
#endif
  !! Display parallel settings
  IF(VERBOSE_LVL .GE. 1) THEN
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    DO i_ = 0,num_procs-1
      CALL mpi_barrier(MPI_COMM_WORLD, ierr)
      IF (my_id .EQ. i_) THEN
        IF (my_id .EQ. 0) WRITE(*,*) '------------ Parallel environment ------------'
        IF (my_id .EQ. 0) WRITE(*,'(A,I4)') ' Total number of processes: ', num_procs
        IF (my_id .EQ. 0) WRITE(*,'(A,I4,A,I4,A,I4,A)') ' Process distribution:'
        IF (my_id .EQ. 0) WRITE(*,'(A,I4,A,I4,A,I4,A)') ' (', num_procs_p, ',', num_procs_ky, ',', num_procs_z,')'
        IF (my_id .EQ. 0) WRITE(*,*) ''
        WRITE(*,'(A10,I4,A3,I4,A1,I4,A1,I4,A1)')&
        ' Process #', my_id, ' (', rank_p, ',', rank_ky,',', rank_z,')'
        WRITE(*,'(A18,I4,A11,I4)')&
        '     local_np   = ', local_np  , ', offset = ', local_np_offset
        WRITE(*,'(A18,I4,A11,I4)')&
        '     local_nj   = ', local_nj  , ', offset = ', local_nj_offset
        WRITE(*,'(A18,I4,A11,I4)')&
        '     local_nkx  = ', local_nkx , ', offset = ', local_nkx_offset
        WRITE(*,'(A18,I4,A11,I4)')&
        '     local_nky  = ', local_nky , ', offset = ', local_nky_offset
        WRITE(*,'(A18,I4,A11,I4)')&
        '     local_nz   = ', local_nz  , ', offset = ', local_nz_offset
        IF (my_id .NE. num_procs-1) WRITE (*,*) ''
        IF (my_id .EQ. num_procs-1) WRITE(*,*) '----------------------------------------------'
        IF (my_id .EQ. num_procs-1) CALL FLUSH(OUTPUT_UNIT)
      ENDIF
    ENDDO
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  ENDIF
  SELECT CASE(hierarchy_closure)
  CASE('truncation')
    CALL speak('Truncation closure',2)
  CASE('max_degree')
    CALL speak('Max degree closure -> Maximal Napj degree is D = '// str(dmax),2)
  END SELECT

END SUBROUTINE auxval
