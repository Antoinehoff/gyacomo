function [s,x, Fs, Fx] = compute_fa_1D(data, options)
%% Compute the dispersion function from the moment hierarchi decomp.
% Normalized Hermite
Hp = @(p,s) polyval(HermitePoly(p),s)./sqrt(2.^p.*factorial(p));
% Hp = @(p,s) hermiteH(p,s)./sqrt(2.^p.*factorial(p));
% Laguerre
Lj = @(j,x) polyval(LaguerrePoly(j),x);
% Maxwellian factor
FaM = @(s,x) exp(-s.^2-x);

s = options.SPAR;
x = options.XPERP;
smin = min(abs(s));
xmin = min(abs(x));

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

% x = 0
Fs = zeros(data.Nkx,data.Nky,numel(s));
FAM = FaM(s,xmin);
for ip_ = 1:Np
    p_ = parray(ip_);
    HH = Hp(p_,s);
    for ij_ = 1:Nj
        j_ = jarray(ij_);
        LL = Lj(j_,xmin);
        HLF = HH.*LL.*FAM;
        for ikx = 1:data.Nkx
            for iky = 1:data.Nky
                Fs(ikx,iky,:) = squeeze(Fs(ikx,iky,:))' + Napj_(ip_,ij_,ikx,iky)*HLF;
            end
        end
    end
end

% s = 0
Fx = zeros(data.Nkx,data.Nky,numel(x));
FAM = FaM(smin,x);
for ip_ = 1:Np
    p_ = parray(ip_);
    HH = Hp(p_,smin);
    for ij_ = 1:Nj
        j_ = jarray(ij_);
        LL = Lj(j_,x);
        HLF = HH.*LL.*FAM;
        for ikx = 1:data.Nkx
            for iky = 1:data.Nky
                Fx(ikx,iky,:) = squeeze(Fx(ikx,iky,:))' + Napj_(ip_,ij_,ikx,iky)*HLF;
            end
        end
    end
end

Fs = real(Fs.*conj(Fs)); %|f_a|^2
Fs = sqrt(squeeze(sum(sum(Fs,1),2))); %sqrt(<|f_a|>kx,ky)
Fs = Fs./max(max(Fs));
Fx = real(Fx.*conj(Fx)); %|f_a|^2
Fx = sqrt(squeeze(sum(sum(Fx,1),2))); %sqrt(<|f_a|>kx,ky)
Fx = Fx./max(max(Fx));
end