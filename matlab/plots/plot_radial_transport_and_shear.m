function [FIGURE] = plot_radial_transport_and_shear(DATA, TAVG_0, TAVG_1)
%Compute steady radial transport
tend = TAVG_1; tstart = TAVG_0;
[~,its0D] = min(abs(DATA.Ts0D-tstart));
[~,ite0D]   = min(abs(DATA.Ts0D-tend));
SCALE = (2*pi/DATA.Nx/DATA.Ny)^2;
gamma_infty_avg = mean(DATA.PGAMMA_RI(its0D:ite0D))*SCALE;
gamma_infty_std = std (DATA.PGAMMA_RI(its0D:ite0D))*SCALE;

FIGURE.fig = figure; FIGURE.FIGNAME = ['ZF_transport_drphi','_',DATA.PARAMS]; set(gcf, 'Position',  [100, 100, 1200, 600])
    subplot(311)
%     yyaxis left
        plot(DATA.Ts0D,DATA.PGAMMA_RI*SCALE,'DisplayName','$\langle n_i \partial_y\phi \rangle_y$'); hold on;
        plot(DATA.Ts0D(its0D:ite0D),ones(ite0D-its0D+1,1)*gamma_infty_avg, '-k',...
            'DisplayName',['$\Gamma^{\infty} = $',num2str(gamma_infty_avg),'$\pm$',num2str(gamma_infty_std)]);
        grid on; set(gca,'xticklabel',[]); ylabel('$\Gamma_x$')
        ylim([0,5*abs(gamma_infty_avg)]); xlim([DATA.Ts0D(1),DATA.Ts0D(end)]);
        title(['$\nu_{',DATA.CONAME,'}=$', num2str(DATA.NU), ', $\kappa_N=$',num2str(DATA.K_N),...
        ', $L=',num2str(DATA.L),'$, $N=',num2str(DATA.Nx),'$, $(P,J)=(',num2str(DATA.PMAXI),',',num2str(DATA.JMAXI),')$,',...
        ' $\mu_{hd}=$',num2str(DATA.MU),', $\Gamma^{\infty} \approx $',num2str(gamma_infty_avg)]);
    %% zonal vs nonzonal energies for phi(t) and shear radial profile
    
    % computation
    Ns3D = numel(DATA.Ts3D);
    [KY, KX] = meshgrid(DATA.ky, DATA.kx);
    Ephi_Z           = zeros(1,Ns3D);
    Ephi_NZ_kgt0     = zeros(1,Ns3D);
    Ephi_NZ_kgt1     = zeros(1,Ns3D);
    Ephi_NZ_kgt2     = zeros(1,Ns3D);
    high_k_phi       = zeros(1,Ns3D);
    dxphi           = zeros(DATA.Nx,DATA.Ny,1,Ns3D);
%     dx2phi           = zeros(DATA.Nx,DATA.Ny,1,Ns3D);
    for it = 1:numel(DATA.Ts3D)
        [amp,ikzf] = max(abs((KX~=0).*DATA.PHI(:,1,1,it)));
        Ephi_NZ_kgt0(it) = squeeze(sum(sum(((sqrt(KX.^2+KY.^2)>0.0).*(KX~=0).*(KY~=0).*(KX.^2+KY.^2).*abs(DATA.PHI(:,:,1,it)).^2))));
        Ephi_NZ_kgt1(it) = squeeze(sum(sum(((sqrt(KX.^2+KY.^2)>1.0).*(KX~=0).*(KY~=0).*(KX.^2+KY.^2).*abs(DATA.PHI(:,:,1,it)).^2))));
        Ephi_NZ_kgt2(it) = squeeze(sum(sum(((sqrt(KX.^2+KY.^2)>2.0).*(KX~=0).*(KY~=0).*(KX.^2+KY.^2).*abs(DATA.PHI(:,:,1,it)).^2))));
        Ephi_Z(it)       = squeeze(sum(sum(((KX~=0).*(KY==0).*(KX.^2).*abs(DATA.PHI(:,:,1,it)).^2))));
        dxphi(:,:,1,it)  = real(fftshift(ifft2(-1i*KX.*(DATA.PHI(:,:,1,it)),DATA.Nx,DATA.Ny)));
%         dx2phi(:,:,1,it) = real(fftshift(ifft2(-KX.^2.*(DATA.PHI(:,:,1,it)),DATA.Nx,DATA.Ny)));
    end
    [~,its3D] = min(abs(DATA.Ts3D-tstart));
    [~,ite3D]   = min(abs(DATA.Ts3D-tend));
    % plot
    subplot(312)
    it0 = 1; itend = Ns3D;
    trange = it0:itend;
    pltx = @(x) x;%-x(1);
    plty = @(x) x./max((x(its3D:ite3D)));
%     plty = @(x) x(500:end)./max(squeeze(x(500:end)));
        yyaxis left
        plot(pltx(DATA.Ts3D(trange)),plty(Ephi_Z(trange)),'DisplayName',['Zonal, ',DATA.CONAME]); hold on;
        ylim([-0.1, 1.5]); ylabel('$E_{Z}$')
        yyaxis right
        plot(pltx(DATA.Ts3D(trange)),plty(Ephi_NZ_kgt0(trange)),'DisplayName','$k_p>0$');
        xlim([DATA.Ts3D(it0), DATA.Ts3D(itend)]);
        ylim([-0.1, 1.5]); ylabel('$E_{NZ}$')
        xlabel('$t c_s/R$'); grid on; set(gca,'xticklabel',[]);% xlim([0 500]);
    subplot(313)
        [TY,TX] = meshgrid(DATA.x,DATA.Ts3D);
        pclr = pcolor(TX,TY,squeeze((mean(dxphi(:,:,1,:),2)))'); set(pclr, 'edgecolor','none'); legend('$\langle \partial_x\phi\rangle_y$') %colorbar;
%         pclr = pcolor(TX,TY,squeeze((mean(dx2phi(:,:,1,:),2)))'); set(pclr, 'edgecolor','none'); legend('$\langle \partial_x^2\phi\rangle_y$') %colorbar;
        caxis([-3 3]);
        xlabel('$t c_s/R$'), ylabel('$x/\rho_s$'); colormap(bluewhitered(256))
end