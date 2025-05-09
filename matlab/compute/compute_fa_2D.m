function [SS,XX,FF,FAM] = compute_fa_2D(data, species, s, x, T)
%% Compute the dispersion function from the moment hierarchi decomp.
% Normalized Hermite
% Hp = @(p,s) polyval(HermitePoly(p),s)./sqrt(2.^p.*factorial(p));
Hp = @(p,s) polyval(HermitePhys_norm(p),s);
% Laguerre
Lj = @(j,x) polyval(LaguerrePoly(j),x);
% Maxwellian factor
FaM = @(s,x) exp(-s.^2-x);

%meshgrid for efficient evaluation
[SS, XX] = meshgrid(s, x);

switch species
    case 'e'
        Napj_     = data.Napjz(2,:,:,:,:);

    case 'i'
        Napj_     = data.Napjz(1,:,:,:,:);
end

parray    = double(data.grids.Parray);
jarray    = double(data.grids.Jarray);

options.SHOW_FLUXSURF = 0;
options.SHOW_METRICS  = 0;
options.SHOW_CURVOP   = 0;
[~, geo_arrays] = plot_metric(data,options);
J  = geo_arrays.Jacobian;
Nz = data.grids.Nz;

tmp_ = 0;
for iz = 1:Nz
    tmp_     =  tmp_ + J(iz)*Napj_(:,:,:,iz,:);
end
Napj_ = tmp_/sum(J(iz));

frames = T;
for it = 1:numel(T)
    [~,frames(it)] = min(abs(T(it)-data.Ts3D)); 
end
frames = unique(frames);

Napj_  = mean(Napj_(:,:,:,:,frames),5);

Napj_  = squeeze(Napj_);


Np = numel(parray); Nj = numel(jarray);

FF = zeros(numel(x),numel(s));
FAM = FaM(SS,XX);
for ip_ = 1:Np
    p_ = parray(ip_);
    HH = Hp(p_,SS);
    for ij_ = 1:Nj
        j_  = jarray(ij_);
        LL  = Lj(j_,XX);
        HLF = HH.*LL.*FAM;
        FF  = FF + Napj_(ip_,ij_)*HLF;
   end
end
FF = (FF.*conj(FF)); %|f_a|^2
FF = sqrt(FF); %<|f_a|>kx,ky
FF = FF./max(FF(:));
end