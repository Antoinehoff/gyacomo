SUBROUTINE stepon
   !   Advance one time step, (num_step=4 for Runge Kutta 4 scheme)
   USE advance_field_routine, ONLY: advance_time_level, advance_moments
   USE basic,                 ONLY: nlend, chrono_advf, chrono_pois,&
      chrono_chck, chrono_clos, chrono_ghst,&
      start_chrono, stop_chrono
   USE closure,               ONLY: apply_closure_model
   USE ghosts,                ONLY: update_ghosts_moments, update_ghosts_EM
   use mpi,                   ONLY: MPI_COMM_WORLD
   USE time_integration,      ONLY: ntimelevel
   USE prec_const,            ONLY: xp, sp
#ifdef TEST_SVD
   USE CLA,                  ONLY: test_svd,filter_sv_moments_ky_pj
#endif
   USE ExB_shear_flow, ONLY: Update_ExB_shear_flow
   IMPLICIT NONE

   INTEGER :: num_step, ierr
   LOGICAL :: mlend

   SUBSTEPS:DO num_step=1,ntimelevel ! eg RK4 compute successively k1, k2, k3, k4
!----- TEST !-----
      ! ! Update the ExB shear flow for the next step
      ! ! This call includes :
      ! !  - the ExB shear value (s(ky)) update for the next time step
      ! !  - the kx grid update
      ! !  - the ExB NL correction factor update (exp(+/- ixkySdts))
      ! !  - (optional) the kernel, poisson op. and ampere op update
      ! CALL Update_ExB_shear_flow(num_step)
!-----END TEST !-----
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
      CALL start_chrono(chrono_advf)
        CALL advance_moments
      CALL stop_chrono(chrono_advf)

      ! Closure enforcement of moments
      CALL start_chrono(chrono_clos)
        CALL apply_closure_model
      CALL stop_chrono(chrono_clos)

      ! Exchanges the ghosts values of N_n+1
      CALL start_chrono(chrono_ghst)
        CALL update_ghosts_moments
      CALL stop_chrono(chrono_ghst)

      ! Update electrostatic potential phi_n = phi(N_n+1) and potential vect psi
      CALL start_chrono(chrono_pois)
        CALL solve_EM_fields
      CALL stop_chrono(chrono_pois)
      ! Update EM ghosts
      CALL start_chrono(chrono_ghst)
        CALL update_ghosts_EM
      CALL stop_chrono(chrono_ghst)

      !-  Check before next step
      CALL start_chrono(chrono_chck)
        CALL checkfield_all()
      CALL stop_chrono(chrono_chck)

      !! TEST SINGULAR VALUE DECOMPOSITION
#ifdef TEST_SVD
      ! CALL test_svd
      CALL filter_sv_moments_ky_pj
#endif
      
      IF( nlend ) EXIT ! exit do loop

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      !----- AFTER: All fields are updated for step = n+1

   END DO SUBSTEPS

CONTAINS
!!!! Basic structure to simplify stepon
   SUBROUTINE assemble_RHS
      USE basic,          ONLY: chrono_mrhs, chrono_sapj, chrono_coll, chrono_grad, chrono_napj, start_chrono, stop_chrono
      USE moments_eq_rhs, ONLY: compute_moments_eq_rhs
      USE collision,      ONLY: compute_Capj
      USE nonlinear,      ONLY: compute_Sapj
      USE processing,     ONLY: compute_nadiab_moments, compute_interp_z, compute_gradients_z
      IMPLICIT NONE
      ! compute auxiliary non adiabatic moments array
      CALL start_chrono(chrono_napj)
        CALL compute_nadiab_moments
      CALL stop_chrono(chrono_napj)

      ! compute gradients and interp
      CALL start_chrono(chrono_grad)
        CALL compute_gradients_z
        CALL compute_interp_z
      CALL stop_chrono(chrono_grad)

      ! compute nonlinear term ("if linear" is included inside)
      CALL start_chrono(chrono_sapj)
        CALL compute_Sapj
      CALL stop_chrono(chrono_sapj)

      ! compute collision term ("if coll, if nu >0" is included inside)
      CALL start_chrono(chrono_coll)
        CALL compute_Capj
      CALL stop_chrono(chrono_coll)

      ! compute the moments equation rhs
      CALL start_chrono(chrono_mrhs)
        CALL compute_moments_eq_rhs
      CALL stop_chrono(chrono_mrhs)
   END SUBROUTINE assemble_RHS

   SUBROUTINE checkfield_all ! Check all the fields for inf or nan
      USE utility,ONLY: is_nan, is_inf
      USE fields, ONLY: phi
      USE MPI
      USE model,            ONLY: LINEARITY, FORCE_SYMMETRY
      IMPLICIT NONE
      LOGICAL :: checkf_
      REAL(sp):: sum_
      ! filtering
      IF(LINEARITY .NE. 'linear') CALL anti_aliasing    ! ensure 0 mode for 2/3 rule
      IF(FORCE_SYMMETRY)          CALL enforce_symmetry ! Enforcing symmetry on kx = 0 (looks useless)

      mlend=.FALSE.
      IF(.NOT.nlend) THEN
         sum_    = REAL(SUM(ABS(phi)),sp)
         checkf_ = (is_nan(sum_,'phi') .OR. is_inf(sum_,'phi'))
         mlend   = (mlend .or. checkf_)
         CALL MPI_ALLREDUCE(mlend, nlend, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
      ENDIF
   END SUBROUTINE checkfield_all

   SUBROUTINE anti_aliasing
      USE fields, ONLY: moments
      USE time_integration, ONLY: updatetlevel
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
                        moments(ia,ip,ij,iky,ikx,iz,updatetlevel) =&
                           AA_x(ikx)* AA_y(iky) * moments(ia,ip,ij,iky,ikx,iz,updatetlevel)
                     ENDDO a
                  ENDDO p
               ENDDO j
            ENDDO ky
         ENDDO kx
      ENDDO z
   END SUBROUTINE anti_aliasing

   SUBROUTINE enforce_symmetry ! Force X(k) = X(N-k)* complex conjugate symmetry
      USE fields, ONLY: phi, psi, moments
      USE time_integration, ONLY: updatetlevel
      USE grid,   ONLY: total_nkx, ikx0,iky0, contains_ky0
      IMPLICIT NONE
      INTEGER :: ikx
      IF ( contains_ky0 ) THEN
         ! moments
         ! z: DO iz = 1,local_nz+ngz
         ! j: DO ij=1,local_nj+ngj
         ! p: DO ip=1,local_np+ngp
         ! a: DO ia=1,local_na
         DO ikx=2,total_nkx/2 !symmetry at ky = 0
            moments(:,:,:,iky0,ikx,:,updatetlevel) = &
               CONJG(moments(:,:,:,iky0,total_nkx+2-ikx,:,updatetlevel))
         END DO
         ! must be real at origin and top right
         moments(:,:,:, iky0,ikx0,:,updatetlevel) = &
            REAL(moments(:,:,:, iky0,ikx0,:,updatetlevel),xp)
         ! ENDDO a
         ! ENDDO p
         ! ENDDO j
         ! ENDDO z
         ! Phi
         DO ikx=2,total_nkx/2 !symmetry at ky = 0
            phi(iky0,ikx,:) = phi(iky0,total_nkx+2-ikx,:)
         END DO
         ! must be real at origin
         phi(iky0,ikx0,:) = REAL(phi(iky0,ikx0,:),xp)
         ! Psi
         DO ikx=2,total_nkx/2 !symmetry at ky = 0
            psi(iky0,ikx,:) = psi(iky0,total_nkx+2-ikx,:)
         END DO
         ! must be real at origin
         psi(iky0,ikx0,:) = REAL(psi(iky0,ikx0,:),xp)
      ENDIF
   END SUBROUTINE enforce_symmetry

END SUBROUTINE stepon
