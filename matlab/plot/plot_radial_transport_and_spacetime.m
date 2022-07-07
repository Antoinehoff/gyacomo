function [FIGURE] = plot_radial_transport_and_spacetime(DATA, OPTIONS)
    %Compute steady radial transport
    tend = OPTIONS.TAVG_1; tstart = OPTIONS.TAVG_0;
    [~,its0D] = min(abs(DATA.Ts0D-tstart));
    [~,ite0D]   = min(abs(DATA.Ts0D-tend));
    [~,its3D] = min(abs(DATA.Ts3D-tstart));
    [~,ite3D]   = min(abs(DATA.Ts3D-tend));
    SCALE = DATA.scale;%(1/DATA.Nx/DATA.Ny)^2;
    Gx_infty_avg = mean(DATA.PGAMMA_RI(its0D:ite0D))*SCALE;
    Gx_infty_std = std (DATA.PGAMMA_RI(its0D:ite0D))*SCALE;
    Qx_infty_avg = mean(DATA.HFLUX_X(its0D:ite0D))*SCALE;
    Qx_infty_std = std (DATA.HFLUX_X(its0D:ite0D))*SCALE;
    disp(['G_x=',sprintf('%2.2f',Gx_infty_avg),'+-',sprintf('%2.2f',Gx_infty_std)]);
    disp(['Q_x=',sprintf('%2.2f',Qx_infty_avg),'+-',sprintf('%2.2f',Qx_infty_std)]);
    f_avg_z      = squeeze(mean(DATA.PHI(:,:,:,:),3));
    [~,ikzf] = max(squeeze(mean(abs(f_avg_z(1,:,its3D:ite3D)),3)));
    ikzf = min([ikzf,DATA.Nky]);
    Ns3D = numel(DATA.Ts3D);
    [KX, KY] = meshgrid(DATA.kx, DATA.ky);
    
    %% computations

    % Compute Gamma from ifft matlab
    Gx = zeros(DATA.Ny,DATA.Nx,numel(DATA.Ts3D));
%     for it = 1:numel(DATA.Ts3D)
%         for iz = 1:DATA.Nz
%             Gx(:,:,it)  = Gx(:,:,it) + ifourier_GENE(-1i*KY.*(DATA.PHI(:,:,iz,it)))...
%                           .*ifourier_GENE(DATA.DENS_I(:,:,iz,it));
%         end
%         Gx(:,:,it)  = Gx(:,:,it)/DATA.Nz;
%     end
    Gx_t_mtlb = squeeze(mean(mean(Gx,1),2)); 
    % Compute Heat flux from ifft matlab
    Qx = zeros(DATA.Ny,DATA.Nx,numel(DATA.Ts3D));
%     for it = 1:numel(DATA.Ts3D)
%         for iz = 1:DATA.Nz
%             Qx(:,:,it)  = Qx(:,:,it) + ifourier_GENE(-1i*KY.*(DATA.PHI(:,:,iz,it)))...
%                           .*ifourier_GENE(DATA.TEMP_I(:,:,iz,it));
%         end
%         Qx(:,:,it)  = Qx(:,:,it)/DATA.Nz;
%     end
    Qx_t_mtlb = squeeze(mean(mean(Qx,1),2)); 
    % zonal vs nonzonal energies for phi(t)

    E_Zmode_SK       = zeros(1,Ns3D);
    E_NZmode_SK      = zeros(1,Ns3D);
    for it = 1:numel(DATA.Ts3D)
        E_Zmode_SK(it)   = squeeze(DATA.ky(ikzf).^2.*abs(squeeze(f_avg_z(ikzf,1,it))).^2);
        E_NZmode_SK(it)  = squeeze(sum(sum(((1+KX.^2+KY.^2).*abs(squeeze(f_avg_z(:,:,it))).^2.*(KY~=0)))));
    end


%% Figure    
mvm = @(x) movmean(x,OPTIONS.NMVA);
    FIGURE.fig = figure; FIGURE.FIGNAME = ['ZF_transport_drphi','_',DATA.PARAMS]; set(gcf, 'Position',  [100, 100, 1000, 600])
    subplot(311)
%     yyaxis left
        plot(mvm(DATA.Ts0D),mvm(DATA.PGAMMA_RI*SCALE),'DisplayName','$\langle n_i \partial_y\phi \rangle_y$'); hold on;
%         plot(mvm(DATA.Ts3D),mvm(Gx_t_mtlb),'DisplayName','matlab comp.'); hold on;
%         plot(DATA.Ts0D(its0D:ite0D),ones(ite0D-its0D+1,1)*Gx_infty_avg, '-k',...
%             'DisplayName',['$\Gamma^{\infty} = $',num2str(Gx_infty_avg),'$\pm$',num2str(Gx_infty_std)]);
%         ylabel('$\Gamma_x$')
%         ylim([0,5*abs(Gx_infty_avg)]); 
%         xlim([DATA.Ts0D(1),DATA.Ts0D(end)]);
%     yyaxis right
        plot(mvm(DATA.Ts0D),mvm(DATA.HFLUX_X*SCALE),'DisplayName','$\langle n_i \partial_y\phi \rangle_y$'); hold on;
%         plot(mvm(DATA.Ts3D),mvm(Qx_t_mtlb),'DisplayName','matlab comp.'); hold on;
        plot(DATA.Ts0D(its0D:ite0D),ones(ite0D-its0D+1,1)*Qx_infty_avg, '-k',...
            'DisplayName',['$Q_x^{\infty} = $',num2str(Qx_infty_avg),'$\pm$',num2str(Qx_infty_std)]);
        ylabel('$Q_x$')  
        ylim([0,5*abs(Qx_infty_avg)]); 
        xlim([DATA.Ts0D(1),DATA.Ts0D(end)]);
    grid on; set(gca,'xticklabel',[]); 
    title({DATA.param_title,...
        ['$\Gamma^{\infty} = $',num2str(Gx_infty_avg),'$, Q^{\infty} = $',num2str(Qx_infty_avg)]});
    %% radial shear radial profile
        % computation
    Ns3D = numel(DATA.Ts3D);
    [KX, KY] = meshgrid(DATA.kx, DATA.ky);
    plt = @(x) mean(x(:,:,:),1);
    kycut = max(DATA.ky);
    kxcut = max(DATA.kx);
    LP = (abs(KY)<kycut).*(abs(KX)<kxcut); %Low pass filter
    
    OPTIONS.NAME = OPTIONS.ST_FIELD;
    OPTIONS.PLAN = 'xy';
    OPTIONS.COMP = 'avg';
    OPTIONS.TIME = DATA.Ts3D;
    OPTIONS.POLARPLOT = 0;
    toplot = process_field(DATA,OPTIONS);
    f2plot = toplot.FIELD;
    dframe = ite3D - its3D;
    clim = max(max(max(abs(plt(f2plot(:,:,:))))));
    subplot(313)
        [TY,TX] = meshgrid(DATA.x,DATA.Ts3D(toplot.FRAMES));
        pclr = pcolor(TX,TY,squeeze(plt(f2plot))'); 
        set(pclr, 'edgecolor','none'); 
        legend(['$\langle ',OPTIONS.ST_FIELD,' \rangle_{y,z}$']) %colorbar;
        caxis(clim*[-1 1]);
        cmap = bluewhitered(256);
        colormap(cmap)
        xlabel('$t c_s/R$'), ylabel('$x/\rho_s$'); 
        if OPTIONS.INTERP
            shading interp
        end
    if strcmp(OPTIONS.ST_FIELD,'Gx')
        subplot(311)
        plot(DATA.Ts3D,squeeze(mean(plt(f2plot),1)));
    end
%% Zonal vs NZonal energies    
    subplot(312)
    it0 = 1; itend = Ns3D;
    trange = toplot.FRAMES;
    plt1 = @(x) x;%-x(1);
    plt2 = @(x) x./max((x(:)));
    toplot = sum(squeeze(plt(f2plot)).^2,1); %ST from before
%     plty = @(x) x(500:end)./max(squeeze(x(500:end)));
        yyaxis left
%         plot(plt1(DATA.Ts3D(trange)),plt2(E_Zmode_SK(trange)),'DisplayName','$k_{zf}^2|\phi_{kzf}|^2$');
        plot(plt1(DATA.Ts3D(trange)),plt2(toplot(:)),'DisplayName','Sum $A^2$');
        ylim([-0.1, 1.5]); ylabel('$E_{Z}$')
        yyaxis right
        plot(plt1(DATA.Ts3D(trange)),plt2(E_NZmode_SK(trange)),'DisplayName','$(1+k^2)|\phi^k|^2$');
        xlim([DATA.Ts3D(it0), DATA.Ts3D(itend)]);
        ylim([-0.1, 1.5]); ylabel('$E_{NZ}$')
        xlabel('$t c_s/R$'); grid on; set(gca,'xticklabel',[]);% xlim([0 500]);
end