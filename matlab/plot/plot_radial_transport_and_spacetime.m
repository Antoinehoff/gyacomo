function [FIGURE] = plot_radial_transport_and_spacetime(DATA, TAVG_0, TAVG_1,stfname,Nmvm)
    %Compute steady radial transport
    tend = TAVG_1; tstart = TAVG_0;
    [~,its0D] = min(abs(DATA.Ts0D-tstart));
    [~,ite0D]   = min(abs(DATA.Ts0D-tend));
    SCALE = (1/DATA.Nx/DATA.Ny)^2;
    gamma_infty_avg = mean(DATA.PGAMMA_RI(its0D:ite0D))*SCALE;
    gamma_infty_std = std (DATA.PGAMMA_RI(its0D:ite0D))*SCALE;
    [~,ikzf] = max(squeeze(mean(abs(squeeze(DATA.PHI(:,1,1,:))),2)));
    Ns3D = numel(DATA.Ts3D);
    [KY, KX] = meshgrid(DATA.ky, DATA.kx);
    
    %% computations

    % Compute Gamma from ifft matlab
    Gx = zeros(DATA.Nx,DATA.Ny,numel(DATA.Ts3D));
    for it = 1:numel(DATA.Ts3D)
        for iz = 1:DATA.Nz
            Gx(:,:,it)  = Gx(:,:,it) + ifourier_GENE(-1i*KY.*(DATA.PHI(:,:,iz,it)),[DATA.Nx,DATA.Ny])...
                          .*ifourier_GENE(DATA.DENS_I(:,:,iz,it),[DATA.Nx,DATA.Ny]);
        end
        Gx(:,:,it)  = Gx(:,:,it)/DATA.Nz;
    end
    Gamma_t_mtlb = squeeze(mean(mean(Gx,1),2)); 
    % zonal vs nonzonal energies for phi(t)

    E_Zmode_SK       = zeros(1,Ns3D);
    E_NZmode_SK      = zeros(1,Ns3D);
    for it = 1:numel(DATA.Ts3D)
        E_Zmode_SK(it)   = squeeze(DATA.ky(ikzf).^2.*abs(DATA.PHI(ikzf,1,1,it))^2);
        E_NZmode_SK(it)  = squeeze(sum(sum(((1+KX.^2+KY.^2).*abs(DATA.PHI(:,:,1,it)).^2.*(KY~=0)))));
    end
    [~,its3D] = min(abs(DATA.Ts3D-tstart));
    [~,ite3D]   = min(abs(DATA.Ts3D-tend));

%% Figure    
mvm = @(x) movmean(x,Nmvm);
    FIGURE.fig = figure; FIGURE.FIGNAME = ['ZF_transport_drphi','_',DATA.PARAMS]; set(gcf, 'Position',  [100, 100, 1200, 600])
    subplot(311)
%     yyaxis left
        plot(mvm(DATA.Ts0D),mvm(DATA.PGAMMA_RI*SCALE),'DisplayName','$\langle n_i \partial_y\phi \rangle_y$'); hold on;
        plot(mvm(DATA.Ts3D),mvm(Gamma_t_mtlb),'DisplayName','matlab comp.'); hold on;
        plot(DATA.Ts0D(its0D:ite0D),ones(ite0D-its0D+1,1)*gamma_infty_avg, '-k',...
            'DisplayName',['$\Gamma^{\infty} = $',num2str(gamma_infty_avg),'$\pm$',num2str(gamma_infty_std)]);
        grid on; set(gca,'xticklabel',[]); ylabel('$\Gamma_x$')
        ylim([0,5*abs(gamma_infty_avg)]); xlim([DATA.Ts0D(1),DATA.Ts0D(end)]);
        title([DATA.param_title,', $\Gamma^{\infty} = $',num2str(gamma_infty_avg),'$\pm$',num2str(gamma_infty_std)]);
    % plot
    subplot(312)
    it0 = 1; itend = Ns3D;
    trange = it0:itend;
    plt1 = @(x) x;%-x(1);
    plt2 = @(x) x./max((x(its3D:ite3D)));
%     plty = @(x) x(500:end)./max(squeeze(x(500:end)));
        yyaxis left
        plot(plt1(DATA.Ts3D(trange)),plt2(E_Zmode_SK(trange)),'DisplayName','$k_{zf}^2|\phi_{kzf}|^2$');
        ylim([-0.1, 1.5]); ylabel('$E_{Z}$')
        yyaxis right
        plot(plt1(DATA.Ts3D(trange)),plt2(E_NZmode_SK(trange)),'DisplayName','$(1+k^2)|\phi^k|^2$');
        xlim([DATA.Ts3D(it0), DATA.Ts3D(itend)]);
        ylim([-0.1, 1.5]); ylabel('$E_{NZ}$')
        xlabel('$t c_s/R$'); grid on; set(gca,'xticklabel',[]);% xlim([0 500]);
    %% radial shear radial profile
        % computation
    Ns3D = numel(DATA.Ts3D);
    [KY, KX] = meshgrid(DATA.ky, DATA.kx);
    plt = @(x) mean(x(:,:,1,:),2);
    kycut = max(DATA.ky);
    kxcut = max(DATA.kx);
    LP = (abs(KY)<kycut).*(abs(KX)<kxcut); %Low pass filter
    switch stfname
        case 'phi'
                phi            = zeros(DATA.Nx,DATA.Ny,1,Ns3D);
                for it = 1:numel(DATA.Ts3D)
                    phi(:,:,1,it)  = ifourier_GENE(DATA.PHI(:,:,1,it),[DATA.Nx,DATA.Ny]);
                end
                f2plot = phi; fname = '$\phi$';
        case 'v_y'
                dxphi            = zeros(DATA.Nx,DATA.Ny,1,Ns3D);
                for it = 1:numel(DATA.Ts3D)
                    dxphi(:,:,1,it)  = ifourier_GENE(-1i*KX.*(DATA.PHI(:,:,1,it)).*LP,[DATA.Nx,DATA.Ny]);
                end
                f2plot = dxphi; fname = '$\langle \partial_x\phi\rangle_y$';
        case 'v_x'
                dyphi            = zeros(DATA.Nx,DATA.Ny,1,Ns3D);
                for it = 1:numel(DATA.Ts3D)
                    dyphi(:,:,1,it)  = ifourier_GENE(1i*KY.*(DATA.PHI(:,:,1,it)).*LP,[DATA.Nx,DATA.Ny]);
                end
                f2plot = dyphi; fname = '$\langle \partial_y\phi\rangle_y$';
        case 'szf'
            dx2phi           = zeros(DATA.Nx,DATA.Ny,1,Ns3D);
            for it = 1:numel(DATA.Ts3D)
                dx2phi(:,:,1,it) = ifourier_GENE(-KX.^2.*(DATA.PHI(:,:,1,it)).*LP,[DATA.Nx,DATA.Ny]);
            end
            f2plot = dx2phi; fname = '$\langle \partial_x^2\phi\rangle_y$';
    end
    clim = max(max(abs(plt(f2plot(:,:,1,its3D:ite3D)))));
    subplot(313)
        [TY,TX] = meshgrid(DATA.x,DATA.Ts3D);
        pclr = pcolor(TX,TY,squeeze(plt(f2plot))'); 
        set(pclr, 'edgecolor','none'); 
        legend(fname) %colorbar;
        caxis(clim*[-1 1]);
        cmap = bluewhitered(256);
        colormap(cmap)
        xlabel('$t c_s/R$'), ylabel('$x/\rho_s$'); 
    if strcmp(stfname,'Gx')
        subplot(311)
        plot(DATA.Ts3D,squeeze(mean(plt(f2plot),1)));
    end
end