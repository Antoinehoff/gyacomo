SUBROUTINE stepon
  !   Advance one time step, (num_step=4 for Runge Kutta 4 scheme)
  USE advance_field_routine, ONLY: advance_time_level, advance_moments
  USE basic,                 ONLY: nlend
  USE closure,               ONLY: apply_closure_model
  USE ghosts,                ONLY: update_ghosts_moments, update_ghosts_EM
  use mpi,                   ONLY: MPI_COMM_WORLD
  USE time_integration,      ONLY: ntimelevel
  USE prec_const,            ONLY: dp
  IMPLICIT NONE

  INTEGER :: num_step, ierr
  LOGICAL :: mlend

   DO num_step=1,ntimelevel ! eg RK4 compute successively k1, k2, k3, k4
   !----- BEFORE: All fields+ghosts are updated for step = n
      ! Compute right hand side from current fields
      ! N_rhs(N_n, nadia_n, phi_n, S_n, Tcoll_n)
      CALL assemble_RHS

      ! ---- step n -> n+1 transition
      ! Advance from updatetlevel to updatetlevel+1 (according to num. scheme)
      CALL advance_time_level
      ! ----

      ! Update moments with the hierarchy RHS (step by step)
      ! N_n+1 = N_n + N_rhs(n)
      CALL advance_moments
      ! Closure enforcement of moments
      CALL apply_closure_model
      ! Exchanges the ghosts values of N_n+1
      CALL update_ghosts_moments

      ! Update electrostatic potential phi_n = phi(N_n+1) and potential vect psi
      CALL solve_EM_fields
      CALL update_ghosts_EM

      !-  Check before next step
      CALL checkfield_all()
      IF( nlend ) EXIT ! exit do loop

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
   !----- AFTER: All fields are updated for step = n+1
   END DO

   CONTAINS
     !!!! Basic structure to simplify stepon
     SUBROUTINE assemble_RHS
       USE moments_eq_rhs, ONLY: compute_moments_eq_rhs
       USE collision,      ONLY: compute_Capj
       USE nonlinear,      ONLY: compute_Sapj
       USE processing,     ONLY: compute_nadiab_moments_z_gradients_and_interp
       IMPLICIT NONE
         ! compute auxiliary non adiabatic moments array, gradients and interp
         CALL compute_nadiab_moments_z_gradients_and_interp
         ! compute nonlinear term ("if linear" is included inside)
         CALL compute_Sapj
         ! compute collision term ("if coll, if nu >0" is included inside)
         CALL compute_Capj
         ! compute the moments equation rhs
         CALL compute_moments_eq_rhs
     END SUBROUTINE assemble_RHS

      SUBROUTINE checkfield_all ! Check all the fields for inf or nan
        USE utility,ONLY: checkelem
        USE basic,  ONLY: t0_checkfield, t1_checkfield, tc_checkfield
        USE fields, ONLY: phi, moments
        USE grid,   ONLY: local_na,local_np,local_nj,local_nky,local_nkx,local_nz,&
                          ngp,ngj,ngz
        USE MPI
        USE time_integration, ONLY: updatetlevel
        USE model,            ONLY: LINEARITY
        IMPLICIT NONE
        LOGICAL :: checkf_
        INTEGER :: ia, ip, ij, iky, ikx, iz
        ! Execution time start
        CALL cpu_time(t0_checkfield)

        IF(LINEARITY .NE. 'linear') CALL anti_aliasing   ! ensure 0 mode for 2/3 rule
        IF(LINEARITY .NE. 'linear') CALL enforce_symmetry ! Enforcing symmetry on kx = 0

        mlend=.FALSE.
        IF(.NOT.nlend) THEN
           z: DO iz = 1,local_nz+ngz
           kx:DO ikx= 1,local_nkx
           ky:DO iky=1,local_nky
             checkf_ = checkelem(phi(iky,ikx,iz),' phi')
             mlend= (mlend .or. checkf_)
             j: DO ij=1,local_nj+ngj
             p: DO ip=1,local_np+ngp
             a: DO ia=1,local_na
                 checkf_ = checkelem(moments(ia,ip,ij,iky,ikx,iz,updatetlevel),' moments')
                 mlend   = (mlend .or. checkf_)
             ENDDO a
             ENDDO p
             ENDDO j
           ENDDO ky
           ENDDO kx
           ENDDO z
           CALL MPI_ALLREDUCE(mlend, nlend, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
        ENDIF
        ! Execution time end
        CALL cpu_time(t1_checkfield)
        tc_checkfield = tc_checkfield + (t1_checkfield - t0_checkfield)
      END SUBROUTINE checkfield_all

      SUBROUTINE anti_aliasing
        USE fields, ONLY: moments
        USE grid,   ONLY: local_na,local_np,local_nj,local_nky,local_nkx,local_nz,&
                          ngp,ngj,ngz, AA_x, AA_y
        IMPLICIT NONE
        INTEGER :: ia, ip, ij, iky, ikx, iz
        z: DO iz = 1,local_nz+ngz
        kx:DO ikx= 1,local_nkx
        ky:DO iky=1,local_nky
        j: DO ij=1,local_nj+ngj
        p: DO ip=1,local_np+ngp
        a: DO ia=1,local_na
                  moments(ia,ip,ij,iky,ikx,iz,:) = AA_x(ikx)* AA_y(iky) * moments(ia,ip,ij,iky,ikx,iz,:)
        ENDDO a
        ENDDO p
        ENDDO j
        ENDDO ky
        ENDDO kx
        ENDDO z
      END SUBROUTINE anti_aliasing

      SUBROUTINE enforce_symmetry ! Force X(k) = X(N-k)* complex conjugate symmetry
        USE fields, ONLY: phi, psi, moments
        USE grid,   ONLY: local_na,local_np,local_nj,total_nkx,local_nz,&
                          ngp,ngj,ngz, ikx0,iky0, contains_ky0
        IMPLICIT NONE
        INTEGER :: ia, ip, ij, ikx, iz
        IF ( contains_ky0 ) THEN
          ! moments
          z: DO iz = 1,local_nz+ngz
          j: DO ij=1,local_nj+ngj
          p: DO ip=1,local_np+ngp
          a: DO ia=1,local_na
                DO ikx=2,total_nkx/2 !symmetry at ky = 0
                  moments(ia,ip,ij,iky0,ikx,iz,:) = CONJG(moments(ia,ip,ij,iky0,total_nkx+2-ikx,iz,:))
                END DO
                ! must be real at origin and top right
                moments(ia,ip,ij, iky0,ikx0,iz,:) = REAL(moments(ia,ip,ij, iky0,ikx0,iz,:),dp)
          ENDDO a
          ENDDO p
          ENDDO j
          ENDDO z
          ! Phi
          DO ikx=2,total_nkx/2 !symmetry at ky = 0
            phi(iky0,ikx,:) = phi(iky0,total_nkx+2-ikx,:)
          END DO
          ! must be real at origin
          phi(iky0,ikx0,:) = REAL(phi(iky0,ikx0,:),dp)
          ! Psi
          DO ikx=2,total_nkx/2 !symmetry at ky = 0
            psi(iky0,ikx,:) = psi(iky0,total_nkx+2-ikx,:)
          END DO
          ! must be real at origin
          psi(iky0,ikx0,:) = REAL(psi(iky0,ikx0,:),dp)
        ENDIF
      END SUBROUTINE enforce_symmetry

END SUBROUTINE stepon
