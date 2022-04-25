function [ omega, FIGURE ] = compute_3D_zpinch_growth_rate(DATA, TRANGE, OPTIONS)
FIGURE.FIGNAME = ['growth_rate_kx=0_ky=0_planes',DATA.PARAMS]; 

t   = DATA.Ts3D;
[~,its] = min(abs(t-TRANGE(1)));
[~,ite] = min(abs(t-TRANGE(end)));
nkx = DATA.Nkx; nky = DATA.Nky; nkz = DATA.Nz;
% Remove AA part
if DATA.Nx > 1
    ikxnz = abs(DATA.PHI(:,1,1,1)) > 0;
else
    ikxnz = abs(DATA.PHI(:,1,1,1)) > -1;
end
ikynz = abs(DATA.PHI(1,:,1,1)) > 0;

phi = fft(DATA.PHI(ikxnz,ikynz,:,:),[],3);

omega = zeros(nkx,nky,nkz);

for iz = 1:nkz
    for iy = 1:nky
        for ix = 1:nkx
            omega(ix,iy,iz) = LinearFit_s(t(its:ite),squeeze(abs(phi(ix,iy,iz,its:ite))));
        end
    end
end

%% plot
kx = DATA.kx; ky = DATA.ky; kz = [(0:nkz/2), (-nkz/2+1):-1];
poskx = kx>=0;
posky = ky>=0;
poskz = kz>=0;

kxeq0 = kx==0;
kzeq0 = kz==0;

omega = omega(poskx,posky,poskz);

FIGURE.fig = figure;
nplots = OPTIONS.kxky + OPTIONS.kzkx + OPTIONS.kzky; 
iplot = 1;

if OPTIONS.kxky
[Y_XY,X_XY] = meshgrid(ky(posky),kx(poskx));
subplot(1,nplots,iplot)
    if ~OPTIONS.keq0
        toplot = squeeze(max(real(omega(:,:,:)),[],3));
        pclr= pcolor(X_XY,Y_XY,toplot); set(pclr,'EdgeColor','none');
        xlabel('$k_x$'); ylabel('$k_y$');
        title('$\max(\gamma)_{kz}$');
    else
        toplot = squeeze(real(omega(:,:,kzeq0)));
        pclr= pcolor(X_XY,Y_XY,toplot); set(pclr,'EdgeColor','none');
        xlabel('$k_x$'); ylabel('$k_y$');
        title('$\gamma(k_z=0)$');
    end
    iplot = iplot + 1;
end
if OPTIONS.kzky
[Y_ZY,Z_ZY] = meshgrid(ky(posky),kz(poskz));
subplot(1,nplots,iplot)
    if ~OPTIONS.keq0
        toplot = squeeze(max(real(omega(:,:,:)),[],1));
        pclr= pcolor(Z_ZY,Y_ZY,toplot'); set(pclr,'EdgeColor','none');
        xlabel('$k_x$'); ylabel('$k_y$');
        title('$\max(\gamma)_{kx}$');
    else
        toplot = squeeze(real(omega(kxeq0,:,:)));
        pclr= pcolor(Z_ZY,Y_ZY,toplot'); set(pclr,'EdgeColor','none');
        xlabel('$k_z$'); ylabel('$k_y$');
        title('$\gamma(k_x=0)$');
    end
    iplot = iplot + 1;
end
if OPTIONS.kzkx
[X_ZX,Z_ZX] = meshgrid(kx(poskx),kz(poskz));
subplot(1,nplots,iplot)
    if ~OPTIONS.keq0
        toplot = squeeze(max(real(omega(:,:,:)),[],2));
        pclr= pcolor(Z_ZX,X_ZX,toplot'); set(pclr,'EdgeColor','none');
        xlabel('$k_z$'); ylabel('$k_x$');
        title('$\max(\gamma)_{ky}$');
    else
        toplot = squeeze(real(omega(:,kyeq0,:)));
        pclr= pcolor(Z_ZY,Y_ZY,toplot'); set(pclr,'EdgeColor','none');
        xlabel('$k_z$'); ylabel('$k_x$');
        title('$\gamma(k_y=0)$');
    end
end
shading interp
colormap(bluewhitered);
end
