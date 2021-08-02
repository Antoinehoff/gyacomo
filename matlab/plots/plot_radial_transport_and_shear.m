%Compute steady radial transport
tend = Ts0D(end); tstart = tend - TAVG;
[~,its0D] = min(abs(Ts0D-tstart));
[~,ite0D]   = min(abs(Ts0D-tend));
SCALE = (2*pi/Nx/Ny)^2;
gamma_infty_avg = mean(PGAMMA_RI(its0D:ite0D))*SCALE;
gamma_infty_std = std (PGAMMA_RI(its0D:ite0D))*SCALE;
% Compute steady shearing rate
tend = Ts3D(end); tstart = tend - TAVG;
[~,its2D] = min(abs(Ts3D-tstart));
[~,ite2D]   = min(abs(Ts3D-tend));
shear_infty_avg = mean(mean(shear_maxx_avgy(:,its2D:ite2D),1));
shear_infty_std = std (mean(shear_maxx_avgy(:,its2D:ite2D),1));
% plots
fig = figure; FIGNAME = ['ZF_transport_drphi','_',PARAMS];set(gcf, 'Position',  [100, 100, 1200, 600])
    subplot(311)
%     yyaxis left
        plot(Ts0D,PGAMMA_RI*SCALE,'DisplayName','$\langle n_i d\phi/dz \rangle_z$'); hold on;
        plot(Ts0D(its0D:ite0D),ones(ite0D-its0D+1,1)*gamma_infty_avg, '-k',...
            'DisplayName',['$\Gamma^{\infty} = $',num2str(gamma_infty_avg),'$\pm$',num2str(gamma_infty_std)]);
        grid on; set(gca,'xticklabel',[]); ylabel('$\Gamma_r$')
        ylim([0,5*abs(gamma_infty_avg)]); xlim([Ts0D(1),Ts0D(end)]);
        title(['$\nu_{',CONAME,'}=$', num2str(NU), ', $\eta=$',num2str(ETAB/ETAN),...
        ', $L=',num2str(L),'$, $N=',num2str(Nx),'$, $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
        ' $\mu_{hd}=$',num2str(MU)]);
        %         
    subplot(312)
        clr      = line_colors(1,:);
        lstyle   = line_styles(1);
%         plt = @(x_) mean(x_,1);
        plt = @(x_) x_(1,:);
        plot(Ts3D,plt(shear_maxx_maxy),'DisplayName','$\max_{r,z}(s_\phi)$'); hold on;
        plot(Ts3D,plt(shear_maxx_avgy),'DisplayName','$\max_{r}\langle s_\phi\rangle_z$'); hold on;
        plot(Ts3D,plt(shear_avgx_maxy),'DisplayName','$\max_{z}\langle s_\phi\rangle_r$'); hold on;
        plot(Ts3D,plt(shear_avgx_avgy),'DisplayName','$\langle s_\phi\rangle_{r,z}$'); hold on;
        plot(Ts3D(its2D:ite2D),ones(ite2D-its2D+1,1)*shear_infty_avg, '-k',...
        'DisplayName',['$s^{\infty} = $',num2str(shear_infty_avg),'$\pm$',num2str(shear_infty_std)]);
        ylim([0,shear_infty_avg*5.0]); xlim([Ts0D(1),Ts0D(end)]);
        grid on; ylabel('Shear amp.');set(gca,'xticklabel',[]);% legend('show');
    subplot(313)
        [TY,TX] = meshgrid(x,Ts3D);
        pclr = pcolor(TX,TY,squeeze((mean(dr2phi(:,:,1,:),2)))'); set(pclr, 'edgecolor','none'); legend('Shear ($\langle \partial_r^2\phi\rangle_z$)') %colorbar; 
        caxis(1*shear_infty_avg*[-1 1]); xlabel('$t c_s/R$'), ylabel('$x/\rho_s$'); colormap(bluewhitered(256))
save_figure