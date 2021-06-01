%% Zonal flow spectral analysis
fig = figure; FIGNAME = ['zonal_flow_spectral_analysis_',PARAMS];
tend = Ts0D(end); tstart = tend-TAVG ;
[~,its0D] = min(abs(Ts0D-tstart));
[~,ite0D]   = min(abs(Ts0D-tend));
[~,its2D] = min(abs(Ts2D-tstart));
[~,ite2D]   = min(abs(Ts2D-tend));
set(gcf, 'Position',  [100, 100, 800, 400])
    % Time series analysis (burst period and time frequencies spectrum)
    subplot(121)
    samplerate = Ts0D(2)-Ts0D(1);
    Y = log(PGAMMA_RI(its0D:ite0D)*(2*pi/Nr/Nz)^2);
    [n,~] = size(Y);
    Yy= fft(Y);
    Pot = Yy .* conj(Yy) / n;
    Freq = (samplerate / n * (1:n))';
    Pot(1) = 0;
    nmax = min(20,round(n/2));
    [amax, itmax] = max(Pot);
    plot((0:nmax-1) , Pot(1:nmax)/amax,'DisplayName','$\Gamma_r(\omega)$');hold on;
    plot([itmax-1,itmax-1],[0,1],'--k', 'DisplayName',['$T_{per}\approx',num2str(round(1/Freq(itmax))),'L_\perp/c_s$']);
    legend('show'); grid on; box on; xlabel('Period number'); yticks([]);
    title('$\Gamma_r$ temporal spectrum')
    % Space analysis (spatial period of ZF)
    subplot(122)
    nmax = 20; n = numel(r);
    [TT,NN] = meshgrid(Ts2D(its2D:ite2D),0:n-1);
    Pot = NN;
    for it = 1:ite2D-its2D+1
        Y = mean(real(drphi(:,:,it)),2);
        Yy = fft(Y);
        [n,~] = size(Yy);
        Pot(:,it) = Yy .* conj(Yy) / n;
    end
    [amax, ikmax] = max(mean(Pot,2));
%     pclr = pcolor(NN(1:nmax,:),TT(1:nmax,:),Pot(1:nmax,:)); set(pclr, 'edgecolor','none'); hold on;
    plot(0:nmax,mean(Pot(1:nmax+1,:),2)/amax,'DisplayName','$\langle\partial_r\phi\rangle_z (k_r)$'); hold on;
    plot([ikmax-1,ikmax-1],[0,1],'--k', 'DisplayName',['$L_z=',num2str(2*pi/kr(ikmax)),'\rho_s$']);
    grid on; box on;
    title('ZF spatial spectrum')
    xlabel('radial mode number');  yticks([]); legend('show')
save_figure

%% Shear and phi amp phase space
fig = figure; FIGNAME = ['phi_shear_phase_space_',PARAMS];
set(gcf, 'Position',  [100, 100, 700, 500])
t1 = Ts2D(end); t0 = 0;
[~,its2D] = min(abs(Ts2D-t0)); [~,ite2D] = min(abs(Ts2D-t1));
scatter(phi_maxr_avgz(its2D:ite2D),shear_maxr_avgz(its2D:ite2D),35,Ts2D(its2D:ite2D),'.',...
    'DisplayName',PARAMS); cbar = colorbar;ylabel(cbar,'$t c_s/\rho_s$','Interpreter','LaTeX')
hold on
xlabel('$\langle \phi \rangle_z^r$'); ylabel('$\langle dV_E/dr \rangle_z^r$')
grid on; title('ES pot. vs Shear phase space')
% plot(phi_avgr_maxz(its2D:ite2D),shear_avgr_maxz(its2D:ite2D),'-')
% plot(phi_maxr_maxz(its2D:ite2D),shear_maxr_maxz(its2D:ite2D),'-')
% plot(phi_avgr_avgz(its2D:ite2D),shear_avgr_avgz(its2D:ite2D),'-')
save_figure

if 0
%% density and phi phase space
fig = figure; FIGNAME = ['phi_ni_phase_space_',PARAMS];
set(gcf, 'Position',  [100, 100, 700, 500])
t1 = Ts2D(end); t0 = 0;
[~,its2D] = min(abs(Ts2D-t0)); [~,ite2D] = min(abs(Ts2D-t1));
scatter3(max(mean(ni00(:,:,its2D:ite2D),2),[],1),phi_maxr_avgz(its2D:ite2D),shear_maxr_avgz(its2D:ite2D),35,Ts2D(its2D:ite2D),'.',...
    'DisplayName',PARAMS); cbar = colorbar;ylabel(cbar,'$t c_s/\rho_s$','Interpreter','LaTeX')
hold on
xlabel('$\langle n_i^{00} \rangle_z^r$'); ylabel('$\langle \phi \rangle_z^r$'); zlabel('$\langle dV_E/dr \rangle_z^r$')
grid on; title('ES pot. vs Shear phase space')
% plot(phi_avgr_maxz(its2D:ite2D),shear_avgr_maxz(its2D:ite2D),'-')
% plot(phi_maxr_maxz(its2D:ite2D),shear_maxr_maxz(its2D:ite2D),'-')
% plot(phi_avgr_avgz(its2D:ite2D),shear_avgr_avgz(its2D:ite2D),'-')
% save_figure
end
%% Non zonal quantities
PHI_NZ = PHI;
PHI_NZ(ikmax-1:ikmax+1,:,:) = 0;

phi_nz    = zeros(Nr,Nz,Ns2D);
for it = 1:numel(Ts2D)
    PH_ = PHI_NZ(:,:,it);
    phi_nz (:,:,it)  = real(fftshift(ifft2((PH_),Nr,Nz)));
end
%%
t0    = 1000;
[~, it02D] = min(abs(Ts2D-t0));
[~, it05D] = min(abs(Ts5D-t0));
skip_ = 10; 
DELAY = 0.005*skip_;
FRAMES_2D = it02D:skip_:numel(Ts2D);
if 0
%% Phi non zonal real space
GIFNAME = ['phi_nz',sprintf('_%.2d',JOBNUM),'_',PARAMS];INTERP = 0;
FIELD = real(phi_nz); X = RR; Y = ZZ; T = Ts2D; FRAMES = FRAMES_2D;
FIELDNAME = '$\phi_{NZ}$'; XNAME = '$r/\rho_s$'; YNAME = '$z/\rho_s$';
create_gif
end