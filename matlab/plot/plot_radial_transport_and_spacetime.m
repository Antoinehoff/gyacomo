function [FIGURE] = plot_radial_transport_and_spacetime(DATA, OPTIONS,CODE)
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
    disp(['G_x=',sprintf('%2.2e',Gx_infty_avg),'+-',sprintf('%2.2e',Gx_infty_std)]);
    disp(['Q_x=',sprintf('%2.2e',Qx_infty_avg),'+-',sprintf('%2.2e',Qx_infty_std)]);
%     disp(['Q_x=',sprintf('%2.2e',Qx_infty_avg),'+-',sprintf('%2.2e',Qx_infty_std)]);
    f_avg_z      = squeeze(mean(DATA.PHI(:,:,:,:),3));
    [~,ikzf] = max(squeeze(mean(abs(f_avg_z(1,:,its3D:ite3D)),3)));
    ikzf = min([ikzf,DATA.Nky]);
    Ns3D = numel(DATA.Ts3D);
    [KX, KY] = meshgrid(DATA.kx, DATA.ky);
    %% error estimation
    DT_       = (tend-tstart)/OPTIONS.NCUT;
    Qx_ee     = zeros(1,OPTIONS.NCUT);
    for i = 1:OPTIONS.NCUT
        [~,its_] = min(abs(DATA.Ts0D - (tstart+(i-1)*DT_)));
        [~,ite_] = min(abs(DATA.Ts0D - (tstart+ i   *DT_)));
        Qx_ee(i) = mean(DATA.HFLUX_X(its_:ite_))*SCALE;
    end
    Qx_avg    = mean(Qx_ee);
    Qx_err    =  std(Qx_ee);
%     disp(['Q_avg=',sprintf('%2.2e',Qx_avg),'+-',sprintf('%2.2e',Qx_err)]);
    %% computations

    % Compute zonal and non zonal energies
    E_Zmode_SK       = zeros(1,Ns3D);
    E_NZmode_SK      = zeros(1,Ns3D);
    for it = 1:numel(DATA.Ts3D)
        E_Zmode_SK(it)   = squeeze(DATA.ky(ikzf).^2.*abs(squeeze(f_avg_z(ikzf,1,it))).^2);
        E_NZmode_SK(it)  = squeeze(sum(sum(((1+KX.^2+KY.^2).*abs(squeeze(f_avg_z(:,:,it))).^2.*(KY~=0)))));
    end
    % Compute thermodynamic entropy Eq.(5) Navarro et al. 2012 PoP
    % 1/2 sum_p sum_j Napj^2(k=0) (avg z)
    switch CODE
        case 'GYACOMO'
        Nipjz = sum(sum(sum(sum(conj(DATA.Nipj).*DATA.Nipj))));
        ff = trapz(DATA.z,Nipjz,5);
        E_TE = 0.5*squeeze(ff);
        % Compute electrostatic energy
        E_ES = zeros(size(DATA.Ts5D));
        bi = sqrt(KX.^2+KY.^2)*DATA.sigma_i*sqrt(2*DATA.tau_i); %argument of the kernel
        for it5D = 1:numel(DATA.Ts5D)
            [~,it3D] = min(abs(DATA.Ts3D-DATA.Ts5D(it5D)));
            for in = 1:DATA.Jmaxi
                Knphi = kernel(in-1,bi).*squeeze(trapz(DATA.z,DATA.PHI(:,:,:,it3D),3));
                Ni0n_z= squeeze(trapz(DATA.z,DATA.Nipj(1,in,:,:,:,it5D),5));
                E_ES(it5D) = 0.5*sum(sum(abs(conj(Knphi).*Ni0n_z)));
            end
        end
        otherwise
            E_TE = 0; E_ES =0; DATA.Ts5D =[0 1];
    end

%% Figure    
clr_ = lines(20);
mvm = @(x) movmean(x,OPTIONS.NMVA);
    FIGURE.fig = figure; FIGURE.FIGNAME = ['ZF_transport_drphi','_',DATA.PARAMS]; %set(gcf, 'Position',  [500, 1000, 1000, 600])
    FIGURE.ax1 = subplot(3,1,1,'parent',FIGURE.fig);
        plot(mvm(DATA.Ts0D),mvm(DATA.PGAMMA_RI*SCALE),'--',...
            'color',clr_((DATA.Pmaxi-1)/2-1,:),...
            'DisplayName',['$\Gamma_x$ ',DATA.paramshort]); hold on;
        plot(mvm(DATA.Ts0D),mvm(DATA.HFLUX_X*SCALE),'-',...
            'color',clr_((DATA.Pmaxi-1)/2-1,:),...
            'DisplayName',['$Q_x$ ',DATA.paramshort]); hold on;
        ylabel('Transport')  
        if(~isnan(Qx_infty_avg))
        plot(DATA.Ts0D(its0D:ite0D),ones(ite0D-its0D+1,1)*Qx_infty_avg, '-k',...
            'DisplayName',['$Q_{avg}=',sprintf('%2.2f',Qx_avg),'\pm',sprintf('%2.2f',Qx_err),'$']); legend('show');
            ylim([0,5*abs(Qx_infty_avg)]); 
        else
        plot(DATA.Ts0D(its0D:ite0D),ones(ite0D-its0D+1,1)*Gx_infty_avg, '-k',...
            'DisplayName',['$\Gamma_{avg}=',sprintf('%2.2f',Gx_infty_avg),'\pm',sprintf('%2.2f',Gx_infty_std),'$']); legend('show');
            try ylim([0,5*abs(Gx_infty_avg)]); 
            catch
            end
        end
        xlim([DATA.Ts0D(1),DATA.Ts0D(end)]);
    grid on; set(gca,'xticklabel',[]); 
    title({DATA.param_title,...
        ['$\Gamma^{\infty} = $',num2str(Gx_infty_avg),'$, Q^{\infty} = $',num2str(Qx_infty_avg)]});
    
 %% Free energy    
    FIGURE.ax2 = subplot(3,1,2,'parent',FIGURE.fig);
    yyaxis left
        plot(DATA.Ts5D,E_TE,'DisplayName','$\epsilon_f$'); hold on;
        ylabel('Entropy');%('$\epsilon_f$')
    yyaxis right
        plot(DATA.Ts5D,E_ES,'DisplayName','$\epsilon_\phi$');
        ylabel('ES energy');%('$\epsilon_\phi$')
        xlim([DATA.Ts5D(1), DATA.Ts5D(end)]);
        xlabel('$t c_s/R$'); grid on; set(gca,'xticklabel',[]);% xlim([0 500]);
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
    FIGURE.ax3 = subplot(3,1,3,'parent',FIGURE.fig);
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
%     subplot(312)
%     it0 = 1; itend = Ns3D;
%     trange = toplot.FRAMES;
%     plt1 = @(x) x;%-x(1);
%     plt2 = @(x) x./max((x(:)));
%     toplot = sum(squeeze(plt(f2plot)).^2,1); %ST from before
% %     plty = @(x) x(500:end)./max(squeeze(x(500:end)));
%         yyaxis left
% %         plot(plt1(DATA.Ts3D(trange)),plt2(E_Zmode_SK(trange)),'DisplayName','$k_{zf}^2|\phi_{kzf}|^2$');
%         plot(plt1(DATA.Ts3D(trange)),plt2(toplot(:)),'DisplayName','Sum $A^2$');
%         ylim([-0.1, 1.5]); ylabel('$E_{Z}$')
%         yyaxis right
%         plot(plt1(DATA.Ts3D(trange)),plt2(E_NZmode_SK(trange)),'DisplayName','$(1+k^2)|\phi^k|^2$');
%         xlim([DATA.Ts3D(it0), DATA.Ts3D(itend)]);
%         ylim([-0.1, 1.5]); ylabel('$E_{NZ}$')
%         xlabel('$t c_s/R$'); grid on; set(gca,'xticklabel',[]);% xlim([0 500]);
end