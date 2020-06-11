
subroutine poisson
  ! Solve poisson equation to get phi

  USE time_integration, ONLY: updatetlevel
  USE array
  USE fields
  USE space_grid

  use prec_const
  IMPLICIT NONE


  INTEGER:: iz

  DO iz=izs,ize
    laplace_rhs(iz) = exp( theta(iz,updatetlevel) ) - 1._dp ! more generally replace 1 by Zeff, "d2phi/dz2 = N-1"
    ! laplace_rhs(iz) = 1.0_dp - exp( theta(iz,updatetlevel) ) ! more generally replace 1 by Zeff, old when using "delec/dz = 1-N"
  END DO
  laplace_rhs(ize+1) = 0._dp ! Boundary condition, insert here desired value of \int_zmin^zmax \phi(z) dz
  ! laplace_rhs(izs+20) = 0.042_dp ! Boundary condition, to fix value of phi(izs) ! OLD, for using phi
  ! laplace_rhs(izs+20) = 0.042_dp ! Boundary condition, to fix mean value of elec
  

!  CALL bsolve(laplace_smumps, laplace_rhs, laplace_sol) ! solves matrix problem "laplace_smumps*laplace_sol=laplace_rhs" for laplace_sol, writes solution in laplace_sol

  DO iz=izs,ize
    phi(iz) =  laplace_sol(iz) ! laplace_sol has an extra coordinate at (ize+1) for the lagrange multiplier
  END DO

END SUBROUTINE poisson
