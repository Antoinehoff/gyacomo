SUBROUTINE moments_eq_rhs  

  USE basic
  USE time_integration
  USE array
  USE fields
  USE fourier_grid
  USE model

  use prec_const
  IMPLICIT NONE

  INTEGER     :: ip,ij, ipj, ipp2j, ipm2j, ipjp1, ipjm1, ikr,ikz
  INTEGER     :: kronp0, kronp2
  INTEGER, DIMENSION(2) :: tmppj
  REAL(dp)    :: ip_dp, ij_dp
  REAL(dp)    :: kz, taueqeetaBi, tauiqietaBi, tauaqaetaBi
  COMPLEX(dp) :: Napj, Napp2j, Napm2j, Napjp1, Napjm1
   
  taueqeetaBi = tau_e/real(q_e)*eta_B * imagu ! Factor tau_s/qu_s*eta_B*i
  tauiqietaBi = tau_i/real(q_i)*eta_B * imagu

  Napp2j = 0; Napm2j = 0;  Napjp1 = 0; Napjm1 = 0; ! higher mixing moments terms
  kronp0 = 0; kronp2 = 0; ! Useful kronecker for phi terms

  ipp2j   = 1; ipm2j   = 1; ipjp1   = 1; ! Mixing moment indices are set to one initially
  ipjm1   = 1;                           ! in order to not overpass the moment array

  do ipj = ipjs, ipje

    if (ipj .le. Nmome + 1) then ! electrons moments

      tmppj   = rabe(ipj) ! Compute p,j electrons moments degree
      ip = tmppj(1); ij = tmppj(2);
      if (ip+2 .le. pmaxe+1) then
        ipp2j   = bare(ip+2,ij)
        Napp2j = moments(ipp2j,ikr,ikz,updatetlevel)
      endif
      if (ip-2 .ge. 1) then
        ipm2j   = bare(ip-2,ij)
        Napm2j = moments(ipm2j,ikr,ikz,updatetlevel)
      endif
      if (ij+1 .le. jmaxe+1) then
        ipjp1   = bare(ip,ij+1)
        Napjp1 = moments(ipjp1,ikr,ikz,updatetlevel)
      endif
      if (ij-1 .ge. 1) then
        ipjm1   = bare(ip,ij-1)
        Napjm1 = moments(ipjm1,ikr,ikz,updatetlevel)
      endif
      tauaqaetaBi = taueqeetaBi !species dependant factor

    else ! ions moments

      tmppj = rabi(ipj) ! Compute p,j ions moments degree
      ip = tmppj(1); ij = tmppj(2);
      if (ip+2 .le. pmaxi+1) then
        ipp2j   = bari(ip+2,ij)
        Napp2j = moments(ipp2j,ikr,ikz,updatetlevel)
      endif
      if (ip-2 .ge. 1) then
        ipm2j   = bari(ip-2,ij)
        Napm2j = moments(ipm2j,ikr,ikz,updatetlevel)
      endif
      if (ij+1 .le. jmaxi+1) then
        ipjp1   = bari(ip,ij+1)
        Napjp1 = moments(ipjp1,ikr,ikz,updatetlevel)
      endif
      if (ij-1 .ge. 1) then
        ipjm1   = bari(ip,ij-1)
        Napjm1 = moments(ipjm1,ikr,ikz,updatetlevel)
      endif
      tauaqaetaBi = tauiqietaBi
    endif

    if (ip .eq. 0) then
      kronp0 = 1
    endif
    if (ip .eq. 2) then
      kronp2 = 1
    endif 

    Napj   = moments(ipj,ikr,ikz,updatetlevel)

    ip_dp = real(ip,dp) 
    ij_dp = real(ij,dp) 

    do ikr = ikrs,ikre 
      do ikz = ikzs,ikze
        kz = kzarray(ikz)
      ! N_a^{pj} term
        Napj   = 2*(ip_dp + ij_dp + 1._dp) * Napj
      ! N_a^{p+2,j} term
        Napp2j = SQRT(ip_dp + 1._dp) * SQRT(ip_dp + 2._dp) * Napp2j
      ! N_a^{p-2,j} term
        Napm2j = SQRT(ip_dp) * SQRT(ip_dp - 1._dp) * Napm2j
      ! N_a^{p,j+1} term
        Napjp1 = -(ij_dp + 1._dp) * Napjp1
      ! N_a^{p,j-1} term
        Napjm1 = -ij_dp * Napjm1

        moments_rhs(ipj,ikr,ikz,updatetlevel) = &
         tauaqaetaBi * kz * (Napj + Napp2j + Napm2j + Napjp1 + Napjm1)
      end do
    end do
  end do

END SUBROUTINE moments_eq_rhs
