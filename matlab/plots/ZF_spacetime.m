function [FIGURE] = ZF_spacetime(DATA, TAVG_0, TAVG_1)
%Compute steady radial transport
tend = TAVG_1; tstart = TAVG_0;
[~,its3D] = min(abs(DATA.Ts3D-tstart));
[~,ite3D]   = min(abs(DATA.Ts3D-tend));
nx = DATA.Nx; ny = DATA.Ny; nz = DATA.Nz; nkx = DATA.Nkx; nky = DATA.Nky;

FIGURE.fig = figure; FIGURE.FIGNAME = ['ZF_spacetime','_',DATA.PARAMS]; set(gcf, 'Position',  [100, 100, 1200, 600])

%% radial shear radial profiles
        % computation
    Ns3D = numel(DATA.Ts3D);
    [KY, KX] = meshgrid(DATA.ky, DATA.kx);
    plt = @(x) mean(x(:,:,1,:),2);
subplot(311)
    phi            = zeros(nx,ny,1,Ns3D);
    for it = 1:numel(DATA.Ts3D)
        % PERSONAL METHOD
%         phi(:,:,1,it)  = real(fftshift(ifft2((DATA.PHI(:,:,1,it)),DATA.Nx,DATA.Ny)));
        % GENE POST PROCESSING METHOD
        phi(:,:,1,it) = ifourier_GENE(DATA.PHI(:,:,1,it),[DATA.Nx,DATA.Ny]);
    end
    clear  spectrumKxKyZ spectrumKxKyZ
    f2plot = phi; fname = '$\langle \phi\rangle_y$';
    [TY,TX] = meshgrid(DATA.x,DATA.Ts3D);
    pclr = pcolor(TX,TY,squeeze(plt(f2plot))'); 
    set(pclr, 'edgecolor','none'); 
    legend(fname) %colorbar;
    clim = max(max(abs(plt(f2plot(:,:,1,its3D:ite3D)))));
    caxis(clim*[-1 1]);
    cmap = bluewhitered(256);
    colormap(cmap); colorbar;
    ylabel('$x/\rho_s$');     
    title(DATA.param_title);
   
subplot(312)
    dxphi            = zeros(DATA.Nx,DATA.Ny,1,Ns3D);
    for it = 1:numel(DATA.Ts3D)
        dxphi(:,:,1,it) = ifourier_GENE(1i*KX.*DATA.PHI(:,:,1,it),[DATA.Nx,DATA.Ny]);
    end
    f2plot = dxphi; fname = '$\langle \partial_x\phi\rangle_y$';
    [TY,TX] = meshgrid(DATA.x,DATA.Ts3D);
    pclr = pcolor(TX,TY,squeeze(plt(f2plot))'); 
    set(pclr, 'edgecolor','none'); 
    legend(fname) %colorbar;
    clim = max(max(abs(plt(f2plot(:,:,1,its3D:ite3D)))));
    caxis(clim*[-1 1]);
    cmap = bluewhitered(256);
    colormap(cmap); colorbar;
    ylabel('$x/\rho_s$');    
    
subplot(313)
    dx2phi           = zeros(DATA.Nx,DATA.Ny,1,Ns3D);
    for it = 1:numel(DATA.Ts3D)
        dx2phi(:,:,1,it) = ifourier_GENE(-KX.^2.*DATA.PHI(:,:,1,it),[DATA.Nx,DATA.Ny]);
    end
    f2plot = dx2phi; 
%     fname = '$\Gamma_x(x)$'; 
    fname = '$\langle \partial_x^2\phi\rangle_y$';
    [TY,TX] = meshgrid(DATA.x,DATA.Ts3D);
    pclr = pcolor(TX,TY,squeeze(plt(f2plot))'); 
    set(pclr, 'edgecolor','none'); 
    legend(fname) %colorbar;
    clim = max(max(abs(plt(f2plot(:,:,1,its3D:ite3D)))));
    caxis(clim*[-1 1]);
    cmap = bluewhitered(256);
    colormap(cmap); colorbar;
    xlabel('$t c_s/R$'), ylabel('$x/\rho_s$'); 
end