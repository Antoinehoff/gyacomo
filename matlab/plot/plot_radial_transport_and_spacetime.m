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
    disp(['Q_avg=',sprintf('%2.2e',Qx_avg),'+-',sprintf('%2.2e',Qx_err)]);
%% Figure    
clr_ = lines(20);
mvm = @(x) movmean(x,OPTIONS.NMVA);
    FIGURE.fig = figure; FIGURE.FIGNAME = ['ZF_transport_drphi','_',DATA.params_string]; %set(gcf, 'Position',  [500, 1000, 1000, 600])
    FIGURE.ax1 = subplot(2,1,1,'parent',FIGURE.fig);
    for ia = 1:DATA.inputs.Na
        plot(mvm(DATA.Ts0D),mvm(DATA.PGAMMA_RI(ia,:)*SCALE),'--',...
            'color',clr_(max(1,(DATA.grids.Np-1)/2+(ia-1)),:),...
            'DisplayName',['$\Gamma_x$ ',DATA.paramshort]); hold on;
        plot(mvm(DATA.Ts0D),mvm(DATA.HFLUX_X(ia,:)*SCALE),'-',...
            'color',clr_(max(1,(DATA.grids.Np-1)/2+(ia-1)),:),...
            'DisplayName',['$Q_x$ ',DATA.paramshort]); hold on;
        ylabel('Transport')  
        if(~isnan(Qx_infty_avg))
        plot(DATA.Ts0D(its0D:ite0D),ones(ite0D-its0D+1,1)*Qx_avg, '-k',...
            'DisplayName',['$Q_{avg}=',sprintf('%2.2f',Qx_avg),'\pm',sprintf('%2.2f',Qx_err),'$']); legend('show');
            ylim([0,5*abs(Qx_infty_avg)]); 
        else
        plot(DATA.Ts0D(its0D:ite0D),ones(ite0D-its0D+1,1)*Gx_infty_avg, '-k',...
            'DisplayName',['$\Gamma_{avg}=',sprintf('%2.2f',Gx_infty_avg),'\pm',sprintf('%2.2f',Gx_infty_std),'$']); legend('show');
            try ylim([0,5*abs(Gx_infty_avg)]); 
            catch
            end
        end
    end
        xlim([DATA.Ts0D(1),DATA.Ts0D(end)]);
    grid on; set(gca,'xticklabel',[]); 
    title({DATA.param_title,...
        ['$\Gamma^{\infty} = $',num2str(Gx_infty_avg),'$, Q^{\infty} = $',num2str(Qx_infty_avg)]});
%% radial shear radial profile
    % computation
    Ns3D = numel(DATA.Ts3D);
    [KX, KY] = meshgrid(DATA.grids.kx, DATA.grids.ky);
    plt = @(x) mean(x(:,:,:),1);
    kycut = max(DATA.grids.ky);
    kxcut = max(DATA.grids.kx);
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
    FIGURE.ax3 = subplot(2,1,2,'parent',FIGURE.fig);
        [TY,TX] = meshgrid(DATA.grids.x,DATA.Ts3D(toplot.FRAMES));
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
    top_title(DATA.paramshort)
end