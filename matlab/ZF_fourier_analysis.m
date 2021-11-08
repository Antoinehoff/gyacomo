%% Zonal flow spectral analysis
fig = figure; FIGNAME = ['zonal_flow_spectral_analysis_',PARAMS];
tend = TAVG_1; tstart = TAVG_0;
[~,its0D] = min(abs(Ts0D-tstart));
[~,ite0D]   = min(abs(Ts0D-tend));
[~,its3D] = min(abs(Ts3D-tstart));
[~,ite3D]   = min(abs(Ts3D-tend));
TAVG = Ts3D(ite3D)-Ts3D(its3D);
set(gcf, 'Position',  [100, 100, 800, 400])
    % Time series analysis (burst period and time frequencies spectrum)
    subplot(121)
    samplerate = Ts0D(2)-Ts0D(1);
    Y = log(PGAMMA_RI(its0D:ite0D)*(2*pi/Nx/Ny)^2);
    [n,~] = size(Y);
    Yy= fft(Y);
    Pot = Yy .* conj(Yy) / n;
    Freq = (samplerate / n * (1:n))';
    Pot(1) = 0;
    nmax = min(20,round(n/2));
    [amax, itmax] = max(Pot);
%     plot((0:nmax-1) , Pot(1:nmax)/amax,'DisplayName','$\Gamma_x(\omega)$');hold on;
%     plot([itmax-1,itmax-1],[0,1],'--k', 'DisplayName',['$T_{per}\approx',num2str(round(1/Freq(itmax))),'L_\perp/c_s$']);
    semilogx(TAVG./(0:nmax-1) , Pot(1:nmax)/amax,'o-','DisplayName','$\Gamma_x(\omega)$');hold on;
    semilogx(TAVG./[itmax-1,itmax-1],[0,1],'--k', 'DisplayName',['$T_{max}\approx',num2str(round(1/Freq(itmax-1))),'L_\perp/c_s$']);
    legend('show'); grid on; box on; xlabel('Period $Tc_s/R$'); yticks([]);
    title('$\Gamma_x$ temporal spectrum')
    % Space analysis (spatial period of ZF)
    subplot(122)
    nmax = 20; n = numel(x);
    [TT,NN] = meshgrid(Ts3D(its3D:ite3D),0:n-1);
    Pot = NN;
    for it = 1:ite3D-its3D+1
        Y = mean(real(dxphi(:,:,it)),2);
        Yy = fft(Y);
        [n,~] = size(Yy);
        Pot(:,it) = Yy .* conj(Yy) / n;
    end
    [amax, ikZF] = max(mean(Pot,2));
%     pclr = pcolor(NN(1:nmax,:),TT(1:nmax,:),Pot(1:nmax,:)); set(pclr, 'edgecolor','none'); hold on;
%     plot(0:nmax,mean(Pot(1:nmax+1,:),2)/amax,'DisplayName','$\langle\partial_x\phi\rangle_y (k_x)$'); hold on;
%     plot([ikZF-1,ikZF-1],[0,1],'--k', 'DisplayName',['$L_x=',num2str(2*pi/kx(ikZF)),'\rho_s$']);
    semilogx(Lx./(0:nmax),mean(Pot(1:nmax+1,:),2)/amax,'o-','DisplayName','$\langle\partial_x\phi\rangle_y (k_x)$'); hold on;
    semilogx(Lx./[ikZF-1,ikZF-1],[0,1],'--k', 'DisplayName',['$L_x=',num2str(2*pi/kx(ikZF)),'\rho_s$']);
    grid on; box on;
    title('ZF spatial spectrum')
    xlabel('Period $\lambda/\rho_s$');  yticks([]); legend('show')
save_figure

%% Pred-Pray phase space (A Zonal Flow review, Diamond 2005, Fig 15, Kobayashi 2015)
% Time evol. of the turbulence energy (Pred in Kobayashi 2015, N = sum phi_k^2 (1+k^2) Non zonal)
E_turb           = zeros(1,Ns3D);
% Time evol. of the ZF energy (Pray in Kobayashi 2015, Ev = phi_q^2 q^2)
E_ZF             = zeros(1,Ns3D);    
for it = 1:numel(Ts3D)
    E_turb(it) = sum(sum(((KY~=0).*(1+KX.^2+KY.^2).*abs(PHI(:,:,1,it)).^2)));
%     E_ZF(it)   = kx(ikZF)^2*abs(PHI(ikZF,1,1,it)).^2;
    E_ZF(it)   = sum(sum(((KY==0).*(1+KX.^2+KY.^2).*abs(PHI(:,:,1,it)).^2)));
end
fig = figure; FIGNAME = ['phi_shear_phase_space_',PARAMS];
set(gcf, 'Position',  [100, 100, 700, 500])
scatter(E_ZF*SCALE,E_turb*SCALE,80,Ts3D,'.',...
    'DisplayName',PARAMS); cbar = colorbar;ylabel(cbar,'$t c_s/\rho_s$','Interpreter','LaTeX')
hold on
% xlabel('$\langle \phi \rangle_z^r$'); ylabel('$\langle dV_E/dr \rangle_z^r$')
xlabel('$E_v$'); ylabel('$N$')
grid on; title('ES pot. vs Shear phase space')
% plot(phi_avgr_maxz(iTs3D:ite2D),shear_avgr_maxz(iTs3D:ite2D),'-')
% plot(phi_maxr_maxz(iTs3D:ite2D),shear_maxr_maxz(iTs3D:ite2D),'-')
% plot(phi_avgr_avgz(iTs3D:ite2D),shear_avgr_avgz(iTs3D:ite2D),'-')
save_figure
clear x_ y_
