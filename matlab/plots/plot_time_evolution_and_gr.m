fig = figure; FIGNAME = ['t_evolutions',sprintf('_%.2d',JOBNUM),'_',PARAMS];
set(gcf, 'Position',  [100, 100, 900, 800])
subplot(111); 
    suptitle(['$\nu_{',CONAME,'}=$', num2str(NU), ', $\eta=$',num2str(ETAB/ETAN),...
        ', $L=',num2str(L),'$, $N=',num2str(Nx),'$, $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
        ' $\mu_{hd}=$',num2str(MU)]);
    subplot(421); 
    for ip = 1:Pe_max
        for ij = 1:Je_max
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$N_e^{',num2str(ip-1),num2str(ij-1),'}$'];
            clr      = line_colors(min(ip,numel(line_colors(:,1))),:);
            lstyle   = line_styles(min(ij,numel(line_styles)));
            semilogy(Ts5D,plt(Ne_norm),'DisplayName',plotname,...
                'Color',clr,'LineStyle',lstyle{1}); hold on;
        end
    end
    grid on; ylabel('$\sum_{k_r,k_z}|N_e^{pj}|$');
    subplot(423)
    for ip = 1:Pi_max
        for ij = 1:Ji_max
            plt      = @(x) squeeze(x(ip,ij,:));
            plotname = ['$N_i^{',num2str(ip-1),num2str(ij-1),'}$'];
            clr      = line_colors(min(ip,numel(line_colors(:,1))),:);
            lstyle   = line_styles(min(ij,numel(line_styles)));
            semilogy(Ts5D,plt(Ni_norm),'DisplayName',plotname,...
                'Color',clr,'LineStyle',lstyle{1}); hold on;
        end
    end
    grid on; ylabel('$\sum_{k_r,k_z}|N_i^{pj}|$'); xlabel('$t c_s/R$')
    subplot(222)
%         plot(Ts0D,GGAMMA_RI*(2*pi/Nx/Ny)^2); hold on;
        yyaxis left; ylabel('$\Gamma_x$');
        plot(Ts3D,squeeze(sum(sum(sum(Gamma_x,1),2),3)));
        yyaxis right; ylabel('$Q_x$');
        plot(Ts3D,squeeze(sum(sum(sum(Q_x,1),2),3)));
        grid on; xlabel('$t c_s/R$');
    if(~isnan(max(max(g_I(1,:,:)))))
    subplot(223)
        plot(ky,max(g_I(1,:,:),[],3),'-','DisplayName','Primar. instability'); hold on;
        plot(kx,max(g_II(:,1,:),[],3),'x-','DisplayName','Second. instability'); hold on;
        plot([max(ky)*2/3,max(ky)*2/3],[0,10],'--k', 'DisplayName','2/3 Orszag AA');
        grid on; xlabel('$k\rho_s$'); ylabel('$\gamma R/c_s$'); legend('show');
        ylim([0,max(max(g_I(1,:,:)))]); xlim([0,max(ky)]);
        shearplot = 426; phiplot = 428;
    else
    shearplot = 223; phiplot = 224;      
    end
    subplot(shearplot)
        plt = @(x) mean(x,1);
        clr      = line_colors(min(ip,numel(line_colors(:,1))),:);
        lstyle   = line_styles(min(ij,numel(line_styles)));
        plot(Ts3D,plt(shear_maxx_maxy),'DisplayName','$\max_{r,z}(s)$'); hold on;
        plot(Ts3D,plt(shear_maxx_avgy),'DisplayName','$\max_{r}\langle s \rangle_z$'); hold on;
        plot(Ts3D,plt(shear_avgx_maxy),'DisplayName','$\max_{z}\langle s \rangle_r$'); hold on;
        plot(Ts3D,plt(shear_avgx_avgy),'DisplayName','$\langle s \rangle_{r,z}$'); hold on;
    grid on; xlabel('$t c_s/R$'); ylabel('$shear$'); 
    subplot(phiplot)
        clr      = line_colors(min(ip,numel(line_colors(:,1))),:);
        lstyle   = line_styles(min(ij,numel(line_styles)));
        plot(Ts3D,plt(phi_maxx_maxy),'DisplayName','$\max_{r,z}(\phi)$'); hold on;
        plot(Ts3D,plt(phi_maxx_avg),'DisplayName','$\max_{r}\langle\phi\rangle_z$'); hold on;
        plot(Ts3D,plt(phi_avgx_maxy),'DisplayName','$\max_{z}\langle\phi\rangle_r$'); hold on;
        plot(Ts3D,plt(phi_avgx_avgy),'DisplayName','$\langle\phi\rangle_{r,z}$'); hold on;
    grid on; xlabel('$t c_s/R$'); ylabel('$E.S. pot$');
save_figure