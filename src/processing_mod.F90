MODULE processing
   USE prec_const,  ONLY: xp, imagu, SQRT2, SQRT3, onetwelfth, twothird
   USE grid,        ONLY: &
      local_na, local_np, local_nj, local_nky, local_nkx, local_nz, Ngz,Ngj,Ngp,nzgrid, &
      parray,pmax,ip0, iodd, ieven,&
      CONTAINSp0,ip1,CONTAINSp1,ip2,CONTAINSp2,ip3,CONTAINSp3,&
      jarray,jmax,ij0, dmax,&
      kyarray, AA_y,&
      kxarray, AA_x,&
      zarray, zweights_SR, ieven, iodd, inv_deltaz
   USE fields,           ONLY: moments, phi, psi
   USE array,            ONLY : kernel, nadiab_moments, &
      ddz_napj, ddzND_Napj, interp_napj,&
      dens, upar, uper, Tpar, Tper, temp
   USE geometry,         ONLY: Jacobian, iInt_Jacobian
   USE time_integration, ONLY: updatetlevel
   USE calculus,         ONLY: simpson_rule_z, grad_z, grad_z_5D, grad_z2, grad_z4, grad_z4_5D, interp_z
   USE model,            ONLY: EM, CLOS, beta, HDz_h
   USE species,          ONLY: tau,q_tau,q_o_sqrt_tau_sigma,sqrt_tau_o_sigma
   USE parallel,         ONLY: num_procs_ky, rank_ky, comm_ky
   USE mpi
   implicit none

   REAL(xp), PUBLIC, ALLOCATABLE, DIMENSION(:), PROTECTED :: pflux_x, gflux_x
   REAL(xp), PUBLIC, ALLOCATABLE, DIMENSION(:), PROTECTED :: hflux_x
   INTEGER :: ierr
   PUBLIC :: init_process
   PUBLIC :: compute_nadiab_moments, compute_gradients_z, compute_interp_z
   PUBLIC :: compute_density, compute_upar, compute_uperp
   PUBLIC :: compute_Tpar, compute_Tperp, compute_fluid_moments
   PUBLIC :: compute_radial_transport
   PUBLIC :: compute_radial_heatflux
   PUBLIC :: compute_Napjz_spectrum
CONTAINS

   SUBROUTINE init_process
      USE grid,       ONLY: local_na
      IMPLICIT NONE
      ALLOCATE( pflux_x(local_na))
      ALLOCATE( gflux_x(local_na))
      ALLOCATE( hflux_x(local_na))
   END SUBROUTINE init_process
!------------------------------ HIGH FREQUENCY ROUTINES -------------
! The following routines (nadiab computing, interp and grads) must be
! as fast as possible since they are called every RK substep.
   ! evaluate the non-adiabatique ion moments
   !
   ! n_{pi} = N^{pj} + kernel(j) /tau_i phi delta_p0
   !
   SUBROUTINE compute_nadiab_moments
      IMPLICIT NONE
      INTEGER :: ia,ip,ij,iky,ikx,iz, j_int, p_int
      !non adiab moments
      DO iz=1,local_nz+ngz
      DO ikx=1,local_nkx
      DO iky=1,local_nky
      DO ij=1,local_nj+ngj
      DO ip=1,local_np+ngp
      DO ia = 1,local_na
         IF(parray(ip) .EQ. 0) THEN
            nadiab_moments(ia,ip,ij,iky,ikx,iz) = moments(ia,ip,ij,iky,ikx,iz,updatetlevel) &
            + q_tau(ia)*kernel(ia,ij,iky,ikx,iz,ieven)*phi(iky,ikx,iz)
            ELSEIF( (parray(ip) .EQ. 1) .AND. (beta .GT. 0) ) THEN
            nadiab_moments(ia,ip,ij,iky,ikx,iz) = moments(ia,ip,ij,iky,ikx,iz,updatetlevel) &
               - q_o_sqrt_tau_sigma(ia)*kernel(ia,ij,iky,ikx,iz,ieven)*psi(iky,ikx,iz)
         ELSE
            nadiab_moments(ia,ip,ij,iky,ikx,iz) = moments(ia,ip,ij,iky,ikx,iz,updatetlevel)
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !! Ensure to kill the moments too high if closue option is set to 1
      IF(CLOS .EQ. 1) THEN
         DO iz=1,local_nz+ngz
         DO ikx=1,local_nkx
         DO iky=1,local_nky
         DO ij=1,local_nj+ngj
         j_int = jarray(ij)
         DO ip=1,local_np+ngp
         p_int = parray(ip)
            DO ia = 1,local_na
            IF(p_int+2*j_int .GT. dmax) &
               nadiab_moments(ia,ip,ij,iky,ikx,iz) = 0._xp
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDIF
   END SUBROUTINE compute_nadiab_moments

   ! z grid gradients
   ! SUBROUTINE compute_gradients_z
   !    IMPLICIT NONE
   !    INTEGER :: eo, p_int, j_int, ia,ip,ij,iky,ikx,iz,izi
   !    COMPLEX(xp), DIMENSION(local_nz+ngz) :: f_in
   !    COMPLEX(xp), DIMENSION(local_nz)     :: f_out
   !       ! Compute z first derivative
   !       DO iz=1,local_nz+ngz
   !          izi = iz+ngz/2
   !       DO ikx=1,local_nkx
   !       DO iky=1,local_nky
   !       DO ij=1,local_nj+ngj
   !       DO ip=1,local_np+ngp
   !       DO ia = 1,local_na
   !          ddz_napj(ia,ip,ij,iky,ikx,iz) = inv_deltaz *(&
   !             +onetwelfth*nadiab_moments(ia,ip,ij,iky,ikx,izi-2)&
   !               -twothird*nadiab_moments(ia,ip,ij,iky,ikx,izi-1)&
   !               +twothird*nadiab_moments(ia,ip,ij,iky,ikx,izi+1)&
   !             -onetwelfth*nadiab_moments(ia,ip,ij,iky,ikx,izi-2)&
   !             )
   !          ddzND_Napj(ia,ip,ij,iky,ikx,iz) = inv_deltaz**4 *(&
   !             +1._xp*moments(ia,ip,ij,iky,ikx,izi-2,updatetlevel)&
   !             -4._xp*moments(ia,ip,ij,iky,ikx,izi-1,updatetlevel)&
   !             +6._xp*moments(ia,ip,ij,iky,ikx,izi  ,updatetlevel)&
   !             -4._xp*moments(ia,ip,ij,iky,ikx,izi+1,updatetlevel)&
   !             +1._xp*moments(ia,ip,ij,iky,ikx,izi-2,updatetlevel)&
   !          )
   !       ENDDO
   !       ENDDO
   !       ENDDO
   !       ENDDO
   !       ENDDO
   !       ENDDO
   ! END SUBROUTINE compute_gradients_z

   ! ! z grid gradients
   SUBROUTINE compute_gradients_z
      IMPLICIT NONE
      INTEGER :: eo, p_int, ia,ip,ij,iky,ikx,iz
      COMPLEX(xp), DIMENSION(local_nz+ngz) :: f_in
      COMPLEX(xp), DIMENSION(local_nz)     :: f_out
      DO ikx = 1,local_nkx
      DO iky = 1,local_nky
      DO ij = 1,local_nj+ngj
      DO ip = 1,local_np+ngp
      DO ia = 1,local_na
         IF(nzgrid .GT. 1) THEN
            p_int = parray(ip+ngp/2)
            eo    = MODULO(p_int,2)+1 ! Indicates if we are on even or odd z grid
         ELSE
            eo = 0
         ENDIF
         ! Compute z first derivative
         f_in = nadiab_moments(ia,ip,ij,iky,ikx,:)
         CALL   grad_z(eo,local_nz,ngz,inv_deltaz,f_in,f_out)
         ddz_napj(ia,ip,ij,iky,ikx,:) = f_out(:)
         ! Parallel numerical diffusion
         IF (HDz_h) THEN
            f_in = nadiab_moments(ia,ip,ij,iky,ikx,:)
         ELSE
            f_in = moments(ia,ip,ij,iky,ikx,:,updatetlevel)
         ENDIF
         CALL  grad_z4(local_nz,ngz,inv_deltaz,f_in,f_out)
         ! get output
         DO iz = 1,local_nz
            ddzND_Napj(ia,ip,ij,iky,ikx,iz) = f_out(iz)
         ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
   END SUBROUTINE compute_gradients_z

      ! z grid interpolation
   SUBROUTINE compute_interp_z
      IMPLICIT NONE
      INTEGER :: eo, ia,ip,ij,iky,ikx,iz
      COMPLEX(xp), DIMENSION(local_nz+ngz) :: f_in
      COMPLEX(xp), DIMENSION(local_nz)     :: f_out
      IF(nzgrid .GT. 1) THEN
         DO ikx = 1,local_nkx
         DO iky = 1,local_nky
         DO ij = 1,local_nj+ngj
         DO ip = 1,local_np+ngp
         DO ia = 1,local_na
            ! Compute even odd grids interpolation
            f_in = nadiab_moments(ia,ip,ij,iky,ikx,:)
            CALL interp_z(eo,local_nz,ngz,f_in,f_out)
            DO iz = 1,local_nz
               interp_napj(ia,ip,ij,iky,ikx,iz) = f_out(iz)
            ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ELSE
         interp_napj(:,:,:,:,:,1:local_nz) = nadiab_moments(:,:,:,:,:,(1+ngz/2):(local_nz+ngz/2))
      ENDIF
   END SUBROUTINE compute_interp_z

   !--------------------- LOW FREQUENCY PROCESSING ROUTINES -------------------!
   ! The following routines are called by the diagnose routine every nsave3D steps 
   ! (does not need to be optimized)
   ! 1D diagnostic to compute the average radial particle transport <n_a v_ExB_x>_xyz
   SUBROUTINE compute_radial_transport
      IMPLICIT NONE
      COMPLEX(xp) :: pflux_local, gflux_local, integral
      REAL(xp)    :: buffer(2)
      INTEGER     :: i_, root, iky, ikx, ia, iz, in, izi, ini
      COMPLEX(xp), DIMENSION(local_nz) :: integrant
      DO ia = 1,local_na
         pflux_local = 0._xp ! particle flux
         gflux_local = 0._xp ! gyrocenter flux
         integrant   = 0._xp ! auxiliary variable for z integration
         !!---------- Gyro center flux (drift kinetic) ------------
         ! Electrostatic part
         IF(CONTAINSp0) THEN
            DO iz = 1,local_nz
            izi = iz + ngz/2 !interior index for ghosted arrays
            DO ikx = 1,local_nkx
            DO iky = 1,local_nky
               integrant(iz) = integrant(iz) &
               +Jacobian(izi,ieven)*moments(ia,ip0,ij0,iky,ikx,izi,updatetlevel)&
               *imagu*kyarray(iky)*CONJG(phi(iky,ikx,izi))
            ENDDO
            ENDDO
            ENDDO
         ENDIF
         ! Electromagnetic part
         IF( EM .AND. CONTAINSp1 ) THEN
            DO iz = 1,local_nz ! we take interior points only
            izi = iz + ngz/2 !interior index for ghosted arrays
            DO ikx = 1,local_nkx
            DO iky = 1,local_nky
               integrant(iz) = integrant(iz)&
                  -Jacobian(izi,iodd)*sqrt_tau_o_sigma(ia)*moments(ia,ip1,ij0,iky,ikx,izi,updatetlevel)&
                  *imagu*kyarray(iky)*CONJG(psi(iky,ikx,izi))
            ENDDO
            ENDDO
            ENDDO
         ENDIF
         ! Integrate over z
         call simpson_rule_z(local_nz,zweights_SR,integrant,integral)
         ! Get process local gyrocenter flux with a factor two to account for the negative ky modes
         gflux_local = 2._xp*integral*iInt_Jacobian
         !
         !!---------- Particle flux (gyrokinetic) ------------
         integrant   = 0._xp ! reset auxiliary variable
         ! Electrostatic part
         IF(CONTAINSp0) THEN
            DO iz = 1,local_nz ! we take interior points only
            izi = iz + ngz/2 !interior index for ghosted arrays
            DO ikx = 1,local_nkx
            DO iky = 1,local_nky
            DO in = 1, local_nj
            ini = in + ngj/2 !interior index for ghosted arrays
               integrant(iz) = integrant(iz)+ &
                  Jacobian(izi,ieven)*moments(ia,ip0,ini,iky,ikx,izi,updatetlevel)&
                  *imagu*kyarray(iky)*kernel(ia,ini,iky,ikx,izi,ieven)*CONJG(phi(iky,ikx,izi))
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDIF
         ! Electromagnetic part
         IF( EM .AND. CONTAINSp1 ) THEN
            DO iz = 1,local_nz ! we take interior points only
            izi = iz + ngz/2 !interior index for ghosted arrays
            DO ikx = 1,local_nkx
            DO iky = 1,local_nky
            DO in = 1, local_nj
            ini = in + ngj/2 !interior index for ghosted arrays
                  integrant(iz) = integrant(iz) - &
                  Jacobian(izi,iodd)*sqrt_tau_o_sigma(ia)*moments(ia,ip1,ini,iky,ikx,izi,updatetlevel)&
                  *imagu*kyarray(iky)*kernel(ia,ini,iky,ikx,izi,iodd)*CONJG(psi(iky,ikx,izi))
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDIF
         ! Integrate over z
         call simpson_rule_z(local_nz,zweights_SR,integrant,integral)
         ! Get process local particle flux with a factor two to account for the negative ky modes
         pflux_local = 2._xp*integral*iInt_Jacobian
         !!!!---------- Sum over all processes ----------
         buffer(1) = REAL(gflux_local,xp)
         buffer(2) = REAL(pflux_local,xp)
         root = 0
         !Gather manually among the rank_p=0 processes and perform the sum
         gflux_x(ia) = 0
         pflux_x(ia) = 0
         IF (num_procs_ky .GT. 1) THEN
            !! Everyone sends its local_sum to root = 0
            IF (rank_ky .NE. root) THEN
               CALL MPI_SEND(buffer, 2 , MPI_DOUBLE_PRECISION, root, 1234, comm_ky, ierr)
            ELSE
               ! Recieve from all the other processes
               DO i_ = 0,num_procs_ky-1
                  IF (i_ .NE. rank_ky) &
                     CALL MPI_RECV(buffer, 2 , MPI_DOUBLE_PRECISION, i_, 1234, comm_ky, MPI_STATUS_IGNORE, ierr)
                  gflux_x(ia) = gflux_x(ia) + buffer(1)
                  pflux_x(ia) = pflux_x(ia) + buffer(2)
               ENDDO
            ENDIF
         ELSE
            gflux_x(ia) = gflux_local
            pflux_x(ia) = pflux_local
         ENDIF
      ENDDO
   END SUBROUTINE compute_radial_transport

! 1D diagnostic to compute the average radial particle heatflux <T_i v_ExB_x>_xyz
   SUBROUTINE compute_radial_heatflux
      IMPLICIT NONE
      COMPLEX(xp) :: hflux_local, integral
      REAL(xp)    :: buffer(2), n_xp
      INTEGER     :: i_, root, in, ia, iky, ikx, iz, izi, ini
      COMPLEX(xp), DIMENSION(local_nz) :: integrant        ! charge density q_a n_a
      DO ia = 1,local_na
         hflux_local = 0._xp ! heat flux
         integrant   = 0._xp ! z integration auxiliary variable
         !!----------------ELECTROSTATIC CONTRIBUTION---------------------------
         IF(CONTAINSp0 .AND. CONTAINSp2) THEN
            ! Loop to compute gamma_kx = sum_ky sum_j -i k_z Kernel_j Na00 * phi
            DO iz = 1,local_nz ! we take interior points only
            izi = iz + ngz/2 !interior index for ghosted arrays
            DO ikx = 1,local_nkx
            DO iky = 1,local_nky
            DO in = 1, local_nj
            ini  = in+ngj/2 !interior index for ghosted arrays
            n_xp = jarray(ini)
               integrant(iz) = integrant(iz) &
                  -Jacobian(izi,ieven)*tau(ia)*imagu*kyarray(iky)*phi(iky,ikx,izi)&
                  *kernel(ia,ini,iky,ikx,izi,ieven)*CONJG(&
                              0.5_xp*SQRT2*moments(ia,ip2,ini  ,iky,ikx,izi,updatetlevel)&
                  +(2._xp*n_xp + 1.5_xp)*moments(ia,ip0,ini  ,iky,ikx,izi,updatetlevel)&
                           -(n_xp+1._xp)*moments(ia,ip0,ini+1,iky,ikx,izi,updatetlevel)&
                                    -n_xp*moments(ia,ip0,ini-1,iky,ikx,izi,updatetlevel))
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ELSEIF(CONTAINSp0) THEN
            ERROR STOP "Degrees p=0 and p=2 should be owned by the same process"
         ENDIF
         IF(EM .AND. CONTAINSp1 .AND. CONTAINSp3) THEN
            !!----------------ELECTROMAGNETIC CONTRIBUTION--------------------
            DO iz  = 1,local_nz
            izi = iz + ngz/2 !interior index for ghosted arrays
            DO ikx = 1,local_nkx
            DO iky = 1,local_nky
            DO in = 1, local_nj
            ini = in + ngj/2 !interior index for ghosted arrays
            n_xp = jarray(ini)
               integrant(iz) = integrant(iz) &
                     +Jacobian(izi,iodd)*tau(ia)*sqrt_tau_o_sigma(ia)*imagu*kyarray(iky)*CONJG(psi(iky,ikx,izi))&
                  *kernel(ia,ini,iky,ikx,izi,iodd)*(&
                  0.5_xp*SQRT2*SQRT3*moments(ia,ip3,ini  ,iky,ikx,izi,updatetlevel)&
                              +1.5_xp*moments(ia,ip1,ini  ,iky,ikx,izi,updatetlevel)&
                  +(2._xp*n_xp+1._xp)*moments(ia,ip1,ini  ,iky,ikx,izi,updatetlevel)&
                        -(n_xp+1._xp)*moments(ia,ip1,ini+1,iky,ikx,izi,updatetlevel)&
                                 -n_xp*moments(ia,ip1,ini-1,iky,ikx,izi,updatetlevel))
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDIF
         ! Add polarisation contribution
         ! integrant(iz) = integrant(iz) + tau_i*imagu*ky_&
         ! *CONJG(phi(iky,ikx,iz))*phi(iky,ikx,iz) * HF_phi_correction_operator(iky,ikx,iz)
         ! Integrate over z
         call simpson_rule_z(local_nz,zweights_SR,integrant,integral)
         ! Double it for taking into account the other half plane
         hflux_local = 2._xp*integral*iInt_Jacobian
         buffer(2)   = REAL(hflux_local,xp)
         root = 0
         !Gather manually among the rank_p=0 processes and perform the sum
         hflux_x(ia) = 0
         IF (num_procs_ky .GT. 1) THEN
            !! Everyone sends its local_sum to root = 0
            IF (rank_ky .NE. root) THEN
               CALL MPI_SEND(buffer, 2 , MPI_DOUBLE_PRECISION, root, 1234, comm_ky, ierr)
            ELSE
               ! Recieve from all the other processes
               DO i_ = 0,num_procs_ky-1
                  IF (i_ .NE. rank_ky) &
                     CALL MPI_RECV(buffer, 2 , MPI_DOUBLE_PRECISION, i_, 1234, comm_ky, MPI_STATUS_IGNORE, ierr)
                  hflux_x(ia) = hflux_x(ia) + buffer(2)
               ENDDO
            ENDIF
         ELSE
            hflux_x(ia) = hflux_local
         ENDIF
      ENDDO
   END SUBROUTINE compute_radial_heatflux

   SUBROUTINE compute_Napjz_spectrum
      USE fields, ONLY : moments
      USE array,  ONLY : Napjz
      USE time_integration, ONLY : updatetlevel
      IMPLICIT NONE
      REAL(xp), DIMENSION(local_np,local_nj,local_nz) :: local_sum,global_sum, buffer
      INTEGER  :: i_, root, count, ia, ip, ij, iky, ikx, iz
      root = 0
      DO ia=1,local_na
         ! z-moment spectrum
         ! build local sum
         local_sum = 0._xp
         DO iz = 1,local_nz
         DO ikx = 1,local_nkx
         DO iky = 1,local_nky
         DO ij = 1,local_nj
         DO ip = 1,local_np
            local_sum(ip,ij,iz)  = local_sum(ip,ij,iz)  + &
               (moments(ia,ip+Ngp/2,ij+Ngj/2,iky,ikx,iz+Ngz/2,updatetlevel) &
               * CONJG(moments(ia,ip+Ngp/2,ij+Ngj/2,iky,ikx,iz+Ngz/2,updatetlevel)))
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ! sum reduction
         buffer     = local_sum
         global_sum = 0._xp
         count = local_np*local_nj*local_nz
         IF (num_procs_ky .GT. 1) THEN
            !! Everyone sends its local_sum to root = 0
            IF (rank_ky .NE. root) THEN
               CALL MPI_SEND(buffer, count , MPI_DOUBLE_PRECISION, root, 5678, comm_ky, ierr)
            ELSE
               ! Recieve from all the other processes
               DO i_ = 0,num_procs_ky-1
                  IF (i_ .NE. rank_ky) &
                     CALL MPI_RECV(buffer, count , MPI_DOUBLE_PRECISION, i_, 5678, comm_ky, MPI_STATUS_IGNORE, ierr)
                  global_sum = global_sum + buffer
               ENDDO
            ENDIF
         ELSE
            global_sum = local_sum
         ENDIF
         Napjz(ia,:,:,:) = global_sum
      ENDDO
   END SUBROUTINE compute_Napjz_spectrum

!_____________________________________________________________________________!
!!!!! FLUID MOMENTS COMPUTATIONS !!!!!
! Compute the 2D particle density for electron and ions (sum over Laguerre)
   SUBROUTINE compute_density
   IMPLICIT NONE
   COMPLEX(xp) :: dens_
   INTEGER :: ia, iz, iky, ikx, ij
   DO ia=1,local_na
   IF ( CONTAINSp0 ) THEN
   ! Loop to compute dens_i = sum_j kernel_j Ni0j at each k
   DO iz = 1,local_nz
   DO iky = 1,local_nky
   DO ikx = 1,local_nkx
      dens_ = 0._xp
      DO ij = 1, local_nj
         dens_ = dens_ + kernel(ia,ij+ngj/2,iky,ikx,iz+ngz/2,ieven) * moments(ia,ip0,ij+ngj/2,iky,ikx,iz+ngz/2,updatetlevel)
      ENDDO
      dens(ia,iky,ikx,iz) = dens_
   ENDDO
   ENDDO
   ENDDO
   ENDIF
   ENDDO
   END SUBROUTINE compute_density

! Compute the 2D particle fluid perp velocity for electron and ions (sum over Laguerre)
   SUBROUTINE compute_uperp
      IMPLICIT NONE
      COMPLEX(xp) :: uperp_
      INTEGER :: ia, iz, iky, ikx, ij
      DO ia=1,local_na
      IF ( CONTAINSp0 ) THEN
      DO iz = 1,local_nz
      DO iky = 1,local_nky
      DO ikx = 1,local_nkx
      uperp_ = 0._xp
      DO ij = 1, local_nj
         uperp_ = uperp_ + kernel(ia,ij+ngj/2,iky,ikx,iz+ngz/2,ieven) *&
            0.5_xp*(moments(ia,ip0,ij+ngj/2,iky,ikx,iz+ngz/2,updatetlevel)&
                     -moments(ia,ip0,ij-1+ngj/2,iky,ikx,iz+ngz/2,updatetlevel))
      ENDDO
      uper(ia,iky,ikx,iz) = uperp_
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      ENDDO
   END SUBROUTINE compute_uperp

! Compute the 2D particle fluid par velocity for electron and ions (sum over Laguerre)
   SUBROUTINE compute_upar
      IMPLICIT NONE
      INTEGER :: ia, iz, iky, ikx, ij
      COMPLEX(xp) :: upar_
      DO ia=1,local_na
      IF ( CONTAINSp1 ) THEN
      DO iz = 1,local_nz
      DO iky = 1,local_nky
      DO ikx = 1,local_nkx
         upar_ = 0._xp
         DO ij = 1, local_nj
            upar_ = upar_ + kernel(ia,ij+ngj/2,iky,ikx,iz+ngz/2,iodd)*moments(ia,ip1,ij+ngj/2,iky,ikx,iz+ngz/2,updatetlevel)
         ENDDO
         upar(ia,iky,ikx,iz) = upar_
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      ENDDO
   END SUBROUTINE compute_upar

! Compute the 2D particle temperature for electron and ions (sum over Laguerre)
   SUBROUTINE compute_tperp
      USE time_integration, ONLY : updatetlevel
      IMPLICIT NONE
      REAL(xp)    :: j_xp
      COMPLEX(xp) :: Tperp_
      INTEGER     :: ia, iz, iky, ikx, ij
      DO ia=1,local_na
      IF ( CONTAINSp0 .AND. CONTAINSp2 ) THEN
      ! Loop to compute T = 1/3*(Tpar + 2Tperp)
      DO iz = 1,local_nz
      DO iky = 1,local_nky
      DO ikx = 1,local_nkx
      Tperp_ = 0._xp
      DO ij = 1, local_nj
         j_xp = jarray(ij+ngj/2)
         Tperp_ = Tperp_ + kernel(ia,ij+ngj/2,iky,ikx,iz+ngz/2,ieven)*&
            ((2_xp*j_xp+1)*moments(ia,ip0,ij  +ngj/2,iky,ikx,iz+ngz/2,updatetlevel)&
                     -j_xp*moments(ia,ip0,ij-1+ngj/2,iky,ikx,iz+ngz/2,updatetlevel)&
                  -(j_xp+1)*moments(ia,ip0,ij+1+ngj/2,iky,ikx,iz+ngz/2,updatetlevel))
      ENDDO
      Tper(ia,iky,ikx,iz) = Tperp_
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      ENDDO
   END SUBROUTINE compute_Tperp

! Compute the 2D particle temperature for electron and ions (sum over Laguerre)
   SUBROUTINE compute_Tpar
      USE time_integration, ONLY : updatetlevel
      IMPLICIT NONE
      REAL(xp)    :: j_xp
      COMPLEX(xp) :: Tpar_
      INTEGER     :: ia, iz, iky, ikx, ij
      DO ia=1,local_na
      IF ( CONTAINSp0 .AND. CONTAINSp0 ) THEN
      ! Loop to compute T = 1/3*(Tpar + 2Tperp)
      DO iz = 1,local_nz
      DO iky = 1,local_nky
      DO ikx = 1,local_nkx
      Tpar_ = 0._xp
      DO ij = 1, local_nj
         j_xp = REAL(ij-1,xp)
         Tpar_  = Tpar_ + kernel(ia,ij+ngj/2,iky,ikx,iz+ngz/2,ieven)*&
            (SQRT2 * moments(ia,ip2,ij+ngj/2,iky,ikx,iz+ngz/2,updatetlevel) &
                     + moments(ia,ip0,ij+ngj/2,iky,ikx,iz+ngz/2,updatetlevel))
      ENDDO
      Tpar(ia,iky,ikx,iz) = Tpar_
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      ENDDO
   END SUBROUTINE compute_Tpar

! Compute the 2D particle fluid moments for electron and ions (sum over Laguerre)
   SUBROUTINE compute_fluid_moments
      IMPLICIT NONE
      CALL compute_density
      CALL compute_upar
      CALL compute_uperp
      CALL compute_Tpar
      CALL compute_Tperp
      ! Temperature
      temp = (Tpar - 2._xp * Tper)/3._xp - dens
   END SUBROUTINE compute_fluid_moments

END MODULE processing
