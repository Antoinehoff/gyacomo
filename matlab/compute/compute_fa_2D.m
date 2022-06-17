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
frames = unique(frames);

Napj_     = mean(Napj_(:,:,:,:,frames),5);

% Napj_ = squeeze(Napj_);


Np = numel(parray); Nj = numel(jarray);

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

if options.RMS
    FF = zeros(data.Nky,data.Nkx,numel(options.XPERP),numel(options.SPAR));
    FAM = FaM(SS,XX);
    for ip_ = 1:Np
        p_ = parray(ip_);
        HH = Hp(p_,SS);
        for ij_ = 1:Nj
            j_ = jarray(ij_);
            LL = Lj(j_,XX);
            HLF = HH.*LL.*FAM;
            for ikx = 1:data.Nkx
                for iky = 1:data.Nky
                    FF(iky,ikx,:,:) = squeeze(FF(iky,ikx,:,:)) + Napj_(ip_,ij_,iky,ikx)*HLF;
                end
            end
       end
    end
else
    FF = zeros(numel(options.XPERP),numel(options.SPAR));
    FAM = FaM(SS,XX);
    for ip_ = 1:Np
        p_ = parray(ip_);
        HH = Hp(p_,SS);
        for ij_ = 1:Nj
            j_ = jarray(ij_);
            LL = Lj(j_,XX);
            FF = FF + squeeze(Napj_(ip_,ij_,ikx0,iky0))*HH.*LL.*FAM;
        end
    end
end
FF = (FF.*conj(FF)); %|f_a|^2
% FF = abs(FF); %|f_a|
if options.RMS
%     FF = squeeze(mean(mean(sqrt(FF),1),2)); %sqrt(<|f_a|^2>kx,ky)
    FF = sqrt(squeeze(mean(mean(FF,1),2))); %<|f_a|>kx,ky
else
    FF = sqrt(squeeze(FF)); %sqrt(<|f_a|>x,y)
end

FF = FF./max(max(FF));
FF = FF';
% FF = sqrt(FF);
% FF = FF';
end