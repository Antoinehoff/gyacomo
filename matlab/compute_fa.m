function FF = compute_fa(Napj, spar, xperp)
%% Compute the dispersion function from the moment hierarchi decomp.
% Normalized Hermite
Hp = @(p,s) polyval(HermitePoly(p),s)./sqrt(2.^p.*factorial(p));
% Laguerre
Lj = @(j,x) polyval(LaguerrePoly(j),x);
% Maxwellian factor
FaM = @(s,x) exp(-s.^2-x);

[SS, XX] = meshgrid(spar, xperp); %meshgrid for efficient evaluation

FF = 0 .* SS;

[Pmax,Jmax] = size(Napj);

FAM = FaM(SS,XX);
for p_ = 0:Pmax-1
    HH = Hp(p_,SS);
    for j_ = 0:Jmax-1
        LL = Lj(j_,XX);
        FF = FF + Napj(p_+1,j_+1)*HH.*LL.*FAM;
    end
end