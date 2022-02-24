function [SS,XX,FF] = compute_fa_2D(data, options)
%% Compute the dispersion function from the moment hierarchi decomp.
% Normalized Hermite
Hp = @(p,s) polyval(HermitePoly(p),s)./sqrt(2.^p.*factorial(p));
% Hp = @(p,s) hermiteH(p,s)./sqrt(2.^p.*factorial(p));
% Laguerre
Lj = @(j,x) polyval(LaguerrePoly(j),x);
% Maxwellian factor
FaM = @(s,x) exp(-s.^2-x);

%meshgrid for efficient evaluation
[SS, XX] = meshgrid(options.SPAR, options.XPERP);

[~, ikx0] = min(abs(data.kx));
[~, iky0] = min(abs(data.ky));

switch options.SPECIE
    case 'e'
        Napj_     = data.Nepj;
        parray    = double(data.Pe);
        jarray    = double(data.Je);
    case 'i'
        Napj_     = data.Nipj;
        parray    = double(data.Pi);
        jarray    = double(data.Ji);
end

switch options.Z
    case 'avg'
        Napj_     = mean(Napj_,5);
    otherwise
        [~,iz]    = min(abs(options.Z-data.z)); 
        Napj_     = Napj_(:,:,:,:,iz,:);
end
Napj_ = squeeze(Napj_);

frames = options.T;
for it = 1:numel(options.T)
    [~,frames(it)] = min(abs(options.T(it)-data.Ts5D)); 
end

Napj_     = mean(Napj_(:,:,:,:,frames),5);

Napj_ = squeeze(Napj_);


Np = numel(parray); Nj = numel(jarray);

FF = zeros(data.Nkx,data.Nky,numel(options.XPERP),numel(options.SPAR));
% FF = zeros(numel(options.XPERP),numel(options.SPAR));

FAM = FaM(SS,XX);
for ip_ = 1:Np
    p_ = parray(ip_);
    HH = Hp(p_,SS);
    for ij_ = 1:Nj
        j_ = jarray(ij_);
        LL = Lj(j_,XX);
%         FF = FF + Napj_(ip_,ij_,ikx0,iky0)*HH.*LL.*FAM;
        HLF = HH.*LL.*FAM;
        for ikx = 1:data.Nkx
            for iky = 1:data.Nky
                FF(ikx,iky,:,:) = squeeze(FF(ikx,iky,:,:)) + Napj_(ip_,ij_,ikx,iky)*HLF;
            end
        end
    end
end
FF = (FF.*conj(FF)); %|f_a|^2
% FF = abs(FF); %|f_a|
FF = sqrt(squeeze(mean(mean(FF,1),2))); %sqrt(<|f_a|>kx,ky)
FF = FF./max(max(FF));
FF = FF';
% FF = FF.*conj(FF);
% FF = sqrt(FF);
% FF = FF./max(max(FF));
% FF = FF';
end