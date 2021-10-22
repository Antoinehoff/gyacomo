%Compute steady radial transport
tend = TAVG_1; tstart = TAVG_0;
[~,its0D] = min(abs(Ts0D-tstart));
[~,ite0D]   = min(abs(Ts0D-tend));
SCALE = (2*pi/Nx/Ny)^2;
gamma_infty_avg = mean(PGAMMA_RI(its0D:ite0D))*SCALE;
gamma_infty_std = std (PGAMMA_RI(its0D:ite0D))*SCALE;
if W_HF
    q_infty_avg     = mean(HFLUX_X(its0D:ite0D))*SCALE;
end
% Compute steady shearing rate
tend = TAVG_1; tstart = TAVG_0;
[~,its2D] = min(abs(Ts3D-tstart));
[~,ite2D]   = min(abs(Ts3D-tend));
shear_infty_avg = mean(mean(shear_maxx_avgy(:,its2D:ite2D),1));
shear_infty_std = std (mean(shear_maxx_avgy(:,its2D:ite2D),1));
% plots
fig = figure; FIGNAME = ['ZF_transport_drphi','_',PARAMS];set(gcf, 'Position',  [100, 100, 1200, 600])
    subplot(311)
    yyaxis left
        plot(Ts0D,PGAMMA_RI*SCALE,'DisplayName','$\langle n_i \partial_y\phi \rangle_y$'); hold on;
        plot(Ts0D(its0D:ite0D),ones(ite0D-its0D+1,1)*gamma_infty_avg, '-k',...
            'DisplayName',['$\Gamma^{\infty} = $',num2str(gamma_infty_avg),'$\pm$',num2str(gamma_infty_std)]);
        grid on; set(gca,'xticklabel',[]); ylabel('$\Gamma_x$')
        ylim([0,5*abs(gamma_infty_avg)]); xlim([Ts0D(1),Ts0D(end)]);
        title(['$\nu_{',CONAME,'}=$', num2str(NU), ', $\kappa_N=$',num2str(K_N),...
        ', $L=',num2str(L),'$, $N=',num2str(Nx),'$, $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
        ' $\mu_{hd}=$',num2str(MU),', $\Gamma^{\infty} \approx $',num2str(gamma_infty_avg)]);
    if W_HF
    yyaxis right
        plot(Ts0D,HFLUX_X*SCALE,'DisplayName','$\langle T_i \partial_y\phi \rangle_y$'); hold on;
        grid on; set(gca,'xticklabel',[]); ylabel('$Q_x$')
        xlim([Ts0D(1),Ts0D(end)]); ylim([0,5*abs(q_infty_avg)]);
    end
        %
%     subplot(312)
%         clr      = line_colors(1,:);
%         lstyle   = line_styles(1);
% %         plt = @(x_) mean(x_,1);
%         plt = @(x_) x_(1,:);
%         plot(Ts3D,plt(shear_maxx_maxy),'DisplayName','$\max_{x,y}(\partial^2_x\phi)$'); hold on;
%         plot(Ts3D,plt(shear_maxx_avgy),'DisplayName','$\max_{x}\langle \partial^2_x\phi\rangle_y$'); hold on;
%         plot(Ts3D,plt(shear_avgx_maxy),'DisplayName','$\max_{y}\langle \partial^2_x\phi\rangle_x$'); hold on;
%         plot(Ts3D,plt(shear_avgx_avgy),'DisplayName','$\langle \partial^2_x\phi\rangle_{x,y}$'); hold on;
%         plot(Ts3D(its2D:ite2D),ones(ite2D-its2D+1,1)*shear_infty_avg, '-k',...
%         'DisplayName',['$s^{\infty} = $',num2str(shear_infty_avg),'$\pm$',num2str(shear_infty_std)]);
%         ylim([0,shear_infty_avg*5.0]); xlim([Ts0D(1),Ts0D(end)]);
%         grid on; ylabel('Shear amp.');set(gca,'xticklabel',[]);% legend('show');
    %% zonal vs nonzonal energies for phi(t)
    subplot(312)
    it0 = 1; itend = Ns3D;
    trange = it0:itend;
    pltx = @(x) x;%-x(1);
    plty = @(x) x./max((x(1:end)));
%     plty = @(x) x(500:end)./max(squeeze(x(500:end)));
        yyaxis left
        plot(pltx(Ts3D(trange)),plty(Ephi_Z(trange)),'DisplayName',['Zonal, ',CONAME]); hold on;
        ylim([1e-1, 1.5]); ylabel('$E_{Z}$')
        yyaxis right
        plot(pltx(Ts3D(trange)),plty(Ephi_NZ_kgt0(trange)),'DisplayName','$k_p>0$');
%         semilogy(pltx(Ts3D(trange)),plty(Ephi_NZ_kgt1(trange)),'DisplayName','$k_p>1$');
%         semilogy(pltx(Ts3D(trange)),plty(Ephi_NZ_kgt2(trange)),'DisplayName','$k_p>2$');
    %     semilogy(pltx(Ts0D),plty(PGAMMA_RI),'DisplayName',['$\Gamma_x$, ',CONAME]);
%         legend('Location','southeast')
        xlim([Ts3D(it0), Ts3D(itend)]);
        ylim([1e-6, 1.0]); ylabel('$E_{NZ}$')
        xlabel('$t c_s/R$'); grid on; set(gca,'xticklabel',[]);% xlim([0 500]);
    subplot(313)
        [TY,TX] = meshgrid(x,Ts3D);
        pclr = pcolor(TX,TY,squeeze((mean(dx2phi(:,:,1,:),2)))'); set(pclr, 'edgecolor','none'); legend('$\langle \partial_x^2\phi\rangle_y$') %colorbar;
        caxis([-3 3]);
        xlabel('$t c_s/R$'), ylabel('$x/\rho_s$'); colormap(bluewhitered(256))
save_figure
