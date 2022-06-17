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

[~, ikx0] = min(abs(data.kx));
[~, iky0] = min(abs(data.ky));
kx_ = data.kx; ky_ = data.ky;

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
Np = numel(parray); Nj = numel(jarray);

switch options.iz
    case 'avg'
        Napj_     = mean(Napj_,5);
        phi_      = mean(data.PHI,3);
    otherwise
        iz        = options.iz; 
        Napj_     = Napj_(:,:,:,:,iz,:);
        phi_      = data.PHI(:,:,iz);
end
% Napj_ = squeeze(Napj_);

frames = options.T;
for it = 1:numel(options.T)
    [~,frames(it)] = min(abs(options.T(it)-data.Ts5D)); 
end

Napj_     = mean(Napj_(:,:,:,:,frames),5);

% Napj_ = squeeze(Napj_);

if options.non_adiab
    for ij_ = 1:Nj
        for ikx = 1:data.Nkx
            for iky = 1:data.Nky    
                kp_ = sqrt(kx_(ikx)^2 + ky_(iky)^2);
                Napj_(1,ij_,iky,ikx) = Napj_(1,ij_,iky,ikx) + kernel(ij_,kp_)*phi_(iky,ikx);
            end
        end
    end
end


% x = 0
if options.RMS
    Fs = zeros(data.Nky,data.Nkx,numel(s));
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
                    Fs(iky,ikx,:) = squeeze(Fs(iky,ikx,:))' + Napj_(ip_,ij_,iky,ikx)*HLF;
                end
            end
       end
    end
else
    Fs = s*0;
    FAM = FaM(s,xmin);
    for ip_ = 1:Np
        p_ = parray(ip_);
        HH = Hp(p_,s);
        for ij_ = 1:Nj
            j_ = jarray(ij_);
            LL = Lj(j_,xmin);
            Fs = Fs + squeeze(Napj_(ip_,ij_,ikx0,iky0))*HH.*LL.*FAM;
        end
    end
end

% s = 0
if options.RMS
    Fx = zeros(data.Nky,data.Nkx,numel(x));
    FAM = FaM(x,smin);
    for ip_ = 1:Np
        p_ = parray(ip_);
        HH = Hp(p_,smin);
        for ij_ = 1:Nj
            j_ = jarray(ij_);
            LL = Lj(j_,x);
            HLF = HH.*LL.*FAM;
            for ikx = 1:data.Nkx
                for iky = 1:data.Nky
                    Fx(iky,ikx,:) = squeeze(Fx(iky,ikx,:))' + Napj_(ip_,ij_,iky,ikx)*HLF;
                end
            end
       end
    end
else
    Fx = x*0;
    FAM = FaM(smin,x);
    for ip_ = 1:Np
        p_ = parray(ip_);
        HH = Hp(p_,smin);
        for ij_ = 1:Nj
            j_ = jarray(ij_);
            LL = Lj(j_,x);
            Fx = Fx + squeeze(Napj_(ip_,ij_,ikx0,iky0))*HH.*LL.*FAM;
        end
    end
end


Fs = real(Fs.*conj(Fs)); %|f_a|^2
Fx = real(Fx.*conj(Fx)); %|f_a|^2
if options.RMS
Fs = squeeze(sqrt(sum(sum(Fs,1),2))); %sqrt(<|f_a|^2>kx,ky)
Fx = squeeze(sqrt(sum(sum(Fx,1),2))); %sqrt(<|f_a|^2>kx,ky)
end
Fs = Fs./max(max(Fs));
Fx = Fx./max(max(Fx));
end