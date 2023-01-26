function [ omega, FIGURE ] = compute_3D_zpinch_growth_rate(DATA, TRANGE, OPTIONS)
FIGURE.FIGNAME = ['growth_rate_kx=0_ky=0_planes',DATA.PARAMS]; 

t   = DATA.Ts3D;
[~,its] = min(abs(t-TRANGE(1)));
[~,ite] = min(abs(t-TRANGE(end)));
nkx = DATA.Nkx; nky = DATA.Nky; nt = numel(t);
% Remove AA part
if DATA.Nx > 1
    ikxnz = abs(DATA.PHI(1,:,1,1)) > 0;
else
    ikxnz = abs(DATA.PHI(1,:,1,1)) > -1;
end
ikynz = abs(DATA.PHI(:,1,1,1)) > 0;

ikynz(1) = 1; ikxnz(1) = 1; %put k=0 in the analysis

% phi = fft(DATA.PHI(:,:,:,:),[],3);
Y = fft(DATA.PHI(ikynz,ikxnz,:,:),[],3);
phi = Y(:,:,2:2:end,:); sz_=size(phi);
nkz = sz_(3);
% tmp_ = ifourier_GENE(DATA.PHI(:,:,:,1)); sz_ = size(tmp_);
% phi  = zeros(sz_(1),sz_(2),sz_(3),nt);
% for it = 1:nt
%    tmp_ = ifourier_GENE(DATA.PHI(:,:,:,it));
%    phi(:,:,:,it) = fftn(tmp_);
% end

omega = zeros(nky,nkx,nkz);

for ikz = 1:nkz
    for iky = 1:nky
        for ikx = 1:nkx
%             omega(iy,ix,iz) = LinearFit_s(t(its:ite),squeeze(abs(phi(iy,ix,iz,its:ite))));
            to_measure = squeeze(log(abs(phi(iky,ikx,ikz,its:ite))));
            tmp = polyfit(t(its:ite),to_measure(:),1);
            if ~(isnan(tmp(1)) || isinf(tmp(1)))
                omega(iky,ikx,ikz) = tmp(1);
            end
        end
    end
end

%% plot
kx = DATA.kx; ky = DATA.ky; kz = [(0:nkz/2), (-nkz/2+1):-1]/DATA.Npol;

kxeq0 = kx==0;
kyeq0 = ky==0;
kzeq0 = kz==0;

% kx = fftshift(kx,1);
% ky = fftshift(ky,1);
% kz = fftshift(kz,1);

FIGURE.fig = figure;
nplots = OPTIONS.kxky + OPTIONS.kzkx + OPTIONS.kzky; 
iplot = 1;

if OPTIONS.kxky
[Y_XY,X_XY] = meshgrid(ky,kx);
subplot(1,nplots,iplot)
    toplot = squeeze(real(omega(:,:,kzeq0)));
    toplot = fftshift(toplot,2);
    X_XY   = fftshift(X_XY,1);
    Y_XY   = fftshift(Y_XY,1);
    pclr= pcolor(X_XY,Y_XY,toplot'); set(pclr,'EdgeColor','none');
    xlabel('$k_x$'); ylabel('$k_y$');
    title('$\gamma(k_z=0)$');
    if OPTIONS.INTERP; shading interp; end
    caxis(max(max(abs(toplot))).*[-1,1]);
    iplot = iplot + 1;
end
if OPTIONS.kzky
[Y_ZY,Z_ZY] = meshgrid(ky,kz);
subplot(1,nplots,iplot)
    toplot = squeeze(real(omega(:,kxeq0,:)));
    toplot = fftshift(toplot,2);
    Z_ZY   = fftshift(Z_ZY,1);
    Y_ZY   = fftshift(Y_ZY,1);
    pclr= pcolor(Z_ZY,Y_ZY,toplot'); set(pclr,'EdgeColor','none');
    xlabel('$k_z$'); ylabel('$k_y$');
    title('$\gamma(k_x=0)$');
    if OPTIONS.INTERP; shading interp; end
    caxis(max(max(abs(toplot))).*[-1,1]);
    iplot = iplot + 1;
end
if OPTIONS.kzkx
[X_ZX,Z_ZX] = meshgrid(kx,kz);
subplot(1,nplots,iplot)
    toplot = squeeze(real(omega(kyeq0,:,:)));
    toplot = fftshift(toplot,2);
    Z_ZX   = fftshift(Z_ZX,1);
    X_ZX   = fftshift(X_ZX,1);
    pclr= pcolor(Z_ZX,X_ZX,toplot'); set(pclr,'EdgeColor','none');
    xlabel('$k_z$'); ylabel('$k_x$');
    title('$\gamma(k_y=0)$');
end
if OPTIONS.INTERP; shading interp; end
caxis(max(max(abs(toplot))).*[-1,1]);
colormap(bluewhitered);
end
