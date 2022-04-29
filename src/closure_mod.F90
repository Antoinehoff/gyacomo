module closure
! Contains the routines to define closures
USE basic
USE model,  ONLY: CLOS, tau_e, tau_i, q_e, q_i, nu, KIN_E
USE grid
USE array,  ONLY: kernel_e,  kernel_i
USE fields, ONLY: moments_e, moments_i
USE time_integration, ONLY: updatetlevel
IMPLICIT NONE

PUBLIC :: apply_closure_model

CONTAINS

! Positive Oob indices are approximated with a model
SUBROUTINE apply_closure_model
  IMPLICIT NONE

  CALL cpu_time(t0_clos)
  IF (CLOS .EQ. 0) THEN
    ! zero truncation, An+1=0 for n+1>nmax only
    CALL ghosts_upper_truncation


  ELSEIF (CLOS .EQ. 1) THEN
    ! Truncation at highest fully represented kinetic moment
    ! e.g. Dmax = 3 means
    ! all Napj s.t. p+2j <= 3
    ! -> (p,j) allowed are (0,0),(1,0),(0,1),(2,0),(1,1),(3,0)
    ! =>> Dmax is Pmax, condition is p+2j<=Pmax
  DO iz = izs,ize
    DO ikx = ikxs,ikxe
        DO iky = ikys,ikye
          IF(KIN_E) THEN
          DO ip = ipgs_e,ipge_e
            DO ij = ijgs_e,ijge_e
              IF ( parray_e(ip)+2*jarray_e(ip) .GT. dmaxe) &
              moments_e(ip,ij,iky,ikx,iz,updatetlevel) = 0._dp
            ENDDO
          ENDDO
          ENDIF
          DO ip = ipgs_i,ipge_i
            DO ij = ijgs_i,ijge_i
              IF ( parray_i(ip)+2*jarray_i(ip) .GT. dmaxi) &
              moments_i(ip,ij,iky,ikx,iz,updatetlevel) = 0._dp
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    ! + ghosts truncation
    CALL ghosts_upper_truncation
  ELSE
    ERROR STOP '! Closure scheme not found !'

  ENDIF

  CALL ghosts_lower_truncation

  CALL cpu_time(t1_clos)
  tc_clos = tc_clos + (t1_clos - t0_clos)
END SUBROUTINE apply_closure_model

! Positive Oob indices are approximated with a model
SUBROUTINE ghosts_upper_truncation
  IMPLICIT NONE

! zero truncation, An+1=0 for n+1>nmax
  ! Electrons
  IF(KIN_E) THEN
    ! applies only for the process that has largest j index
    IF(ije_e .EQ. Jmaxe+1) THEN
      DO ip = ipgs_e,ipge_e
        moments_e(ip,ije_e+1,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel) = 0._dp
      ENDDO
    ENDIF
    ! applies only for the process that has largest p index
    IF(ipe_e .EQ. Pmaxe+1) THEN
      DO ij = ijgs_e,ijge_e
        moments_e(ipe_e+1,ij,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel) = 0._dp
        IF(deltape .EQ. 1) THEN ! Must truncate the second stencil
        moments_e(ipe_e+2,ij,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel) = 0._dp
        ENDIF
      ENDDO
    ENDIF
  ENDIF

  ! Ions
  ! applies only for the process that has largest j index
  IF(ije_i .EQ. Jmaxi+1) THEN
    DO ip = ipgs_i,ipge_i
      moments_i(ip,ije_i+1,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel) = 0._dp
    ENDDO
  ENDIF
  ! applies only for the process that has largest p index
  IF(ipe_i .EQ. Pmaxi+1) THEN
    DO ij = ijgs_i,ijge_i
      moments_i(ipe_i+1,ij,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel) = 0._dp
      IF(deltape .EQ. 1) THEN ! Must truncate the second stencil
      moments_i(ipe_i+2,ij,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel) = 0._dp
      ENDIF
    ENDDO
  ENDIF

END SUBROUTINE ghosts_upper_truncation

! Negative OoB indices are 0
SUBROUTINE ghosts_lower_truncation
  IMPLICIT NONE

! zero truncation, An=0 for n<0
  ! Electrons
  IF(KIN_E) THEN
    ! applies only for the process that has lowest j index
    IF(ijs_e .EQ. 1) THEN
      DO ip = ipgs_e,ipge_e
        moments_e(ip,ijs_e-1,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel) = 0._dp
      ENDDO
    ENDIF
    ! applies only for the process that has lowest p index
    IF(ips_e .EQ. 1) THEN
      DO ij = ijgs_e,ijge_e
        moments_e(ips_e-1,ij,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel) = 0._dp
        IF(deltape .EQ. 1) THEN ! Must truncate the second stencil
        moments_e(ips_e-2,ij,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel) = 0._dp
        ENDIF
      ENDDO
    ENDIF
  ENDIF

  ! Ions
  IF(ijs_i .EQ. 1) THEN
    ! applies only for the process that has lowest j index
    DO ip = ipgs_i,ipge_i
      moments_i(ip,ijs_i-1,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel) = 0._dp
    ENDDO
  ENDIF
  ! applies only for the process that has lowest p index
  IF(ips_i .EQ. 1) THEN
    DO ij = ijgs_i,ijge_i
      moments_i(ips_i-1,ij,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel) = 0._dp
      IF(deltape .EQ. 1) THEN ! Must truncate the second stencil
      moments_i(ips_i-2,ij,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel) = 0._dp
      ENDIF
    ENDDO
  ENDIF

END SUBROUTINE ghosts_lower_truncation

END module closure
