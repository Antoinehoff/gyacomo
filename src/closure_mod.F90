module closure
! Contains the routines to define closures
IMPLICIT NONE

PUBLIC :: apply_closure_model

CONTAINS

! Positive Oob indices are approximated with a model
SUBROUTINE apply_closure_model
  USE basic,      ONLY: t0_clos, t1_clos, tc_clos
  USE prec_const, ONLY: dp
  USE model,      ONLY: CLOS
  USE grid,       ONLY: local_na, local_nky, local_nkx, local_nz,ngz,&
                        local_nj,ngj, jarray,&
                        local_np,ngp, parray, dmax
  USE fields,     ONLY: moments
  USE time_integration, ONLY: updatetlevel
  IMPLICIT NONE
  INTEGER :: iz,ikx,iky,ij,ip,ia
  CALL cpu_time(t0_clos)
  IF (CLOS .EQ. 0) THEN
    ! zero truncation, An+1=0 for n+1>nmax only
    CALL ghosts_upper_truncation
  ELSEIF (CLOS .EQ. 1) THEN
    ! truncation at highest fully represented kinetic moment
    ! e.g. Dmax = 3 means
    ! only Napj s.t. p+2j <= 3 are evolved
    ! -> (p,j) allowed are (0,0),(1,0),(0,1),(2,0),(1,1),(3,0)
    ! =>> Dmax = min(Pmax,2*Jmax+1)
    ! z: DO iz = 1,local_nz+ngz
    ! kx:DO ikx= 1,local_nkx
    ! ky:DO iky= 1,local_nky
    j: DO ij = 1,local_nj+ngj
    p: DO ip = 1,local_np+ngp
      IF ( parray(ip)+2*jarray(ij) .GT. dmax) THEN
        ! a:DO ia = 1,local_na
        ! moments(ia,ip,ij,iky,ikx,iz,updatetlevel) = 0._dp
        moments(ia,ip,ij,:,:,:,updatetlevel) = 0._dp
      ! ENDDO a
      ENDIF
    ENDDO p
    ENDDO j
    ! ENDDO ky
    ! ENDDO kx
    ! ENDDO z
  ELSE
    ERROR STOP '>> ERROR << Closure scheme not found '
  ENDIF
  CALL ghosts_lower_truncation
  CALL cpu_time(t1_clos)
  tc_clos = tc_clos + (t1_clos - t0_clos)
END SUBROUTINE apply_closure_model

! Positive Oob indices are approximated with a model
SUBROUTINE ghosts_upper_truncation
  USE prec_const, ONLY: dp
  USE grid,       ONLY: local_na, local_np,ngp,Pmax,&
                        local_nj,ngj,Jmax,&
                        local_nky,local_nkx,&
                        local_nz,ngz,&
                        local_pmax, local_jmax
  USE fields,           ONLY: moments
  USE time_integration, ONLY: updatetlevel
  IMPLICIT NONE
  INTEGER :: iz,ikx,iky,ip,ij,ia,ig
  ! zero truncation, An+1=0 for n+1>nmax
    ! applies only for the processes that evolve the highest moment degree
    IF(local_jmax .GE. Jmax) THEN
      ! DO iz = 1,local_nz+ngz
      ! DO ikx= 1,local_nkx
      ! DO iky= 1,local_nky
      DO ig = 1,ngj/2
      ! DO ip = 1,local_np+ngp
      ! DO ia = 1,local_na
        ! moments(ia,ip,local_nj+ngj/2+ig,iky,ikx,iz,updatetlevel) = 0._dp
        moments(:,:,local_nj+ngj/2+ig,:,:,:,updatetlevel) = 0._dp
      ! ENDDO
      ! ENDDO
      ENDDO
      ! ENDDO
      ! ENDDO
      ! ENDDO
    ENDIF
    ! applies only for the process that has largest p index
    IF(local_pmax .GE. Pmax) THEN
      ! DO iz  = 1,local_nz+ngz
      ! DO ikx = 1,local_nkx
      ! DO iky = 1,local_nky
      ! DO ij  = 1,local_nj+ngj
      DO ig = 1,ngp/2
        ! DO ia  = 1,local_na
        !   moments(ia,local_np+ngp/2+ig,ij,iky,ikx,iz,updatetlevel) = 0._dp
        moments(:,local_np+ngp/2+ig,:,:,:,:,updatetlevel) = 0._dp
        ! ENDDO
      ENDDO
      ! ENDDO
      ! ENDDO
      ! ENDDO
      ! ENDDO
    ENDIF
END SUBROUTINE ghosts_upper_truncation

! Negative OoB indices are 0
SUBROUTINE ghosts_lower_truncation
  USE prec_const, ONLY: dp
  USE grid,       ONLY: local_na,local_np,ngp,&
                        local_nj,ngj,&
                        local_nky,local_nkx,&
                        local_nz,ngz,&
                        local_pmin, local_jmin
  USE fields,           ONLY: moments
  USE time_integration, ONLY: updatetlevel
  IMPLICIT NONE
  INTEGER :: iz,ikx,iky,ip,ia,ij,ig
! zero truncation, An=0 for n<0
    IF(local_jmin .EQ. 0) THEN
      ! DO iz  = 1,local_nz+ngz
      ! DO ikx = 1,local_nkx
      ! DO iky = 1,local_nky
      DO ig  = 1,ngj/2
      ! DO ip  = 1,local_np+ngp
      ! DO ia  = 1,local_na
        ! set to zero the first ghosts cells
        ! moments(ia,ip,ig,iky,ikx,iz,updatetlevel) = 0._dp
        moments(:,:,ig,:,:,:,updatetlevel) = 0._dp
      ! ENDDO
      ! ENDDO
      ENDDO
      ! ENDDO
      ! ENDDO
      ! ENDDO
    ENDIF
    ! applies only for the process that has lowest p index
    IF(local_pmin .EQ. 0) THEN
      ! DO iz  = 1,local_nz+ngz
      ! DO ikx = 1,local_nkx
      ! DO iky = 1,local_nky
      ! DO ij  = 1,local_nj+ngj
      DO ig  = 1,ngp/2
      ! DO ia  = 1,local_na
        ! moments(ia,ig,ij,iky,ikx,iz,updatetlevel) = 0._dp
        moments(:,ig,:,:,:,:,updatetlevel) = 0._dp
      ENDDO
      ! ENDDO
      ! ENDDO
      ! ENDDO
      ! ENDDO
      ! ENDDO
    ENDIF

END SUBROUTINE ghosts_lower_truncation

END module closure
