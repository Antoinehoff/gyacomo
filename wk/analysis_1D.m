default_plots_options
%% load resulTs2D
JOBNUM = 00;
load_results
Ni00 = squeeze(Ni00);
Ne00 = squeeze(Ne00);
PHI  = squeeze(PHI);
Ts2D      = Ts2D';
Ns      = numel(Ts2D);
dt_samp = mean(diff(Ts2D));
% Build grids
Nkr = numel(kr); Nkz = numel(kz);
[KZ,KR] = meshgrid(kz,kr);
Lkr = max(kr)-min(kr); Lkz = max(kz)-min(kz);
dkr = Lkr/(Nkr-1); dkz = Lkz/(Nkz-1);
Lk = max(Lkr,Lkz);

dr = 2*pi/Lk; dz = 2*pi/Lk;
Nr = max(Nkr,Nkz) * (Nkr > 1) + (Nkr == 1);  Nz = Nkz;
r = dr*(-Nr/2:(Nr/2-1)); Lr = max(r)-min(r);
z = dz*(-Nz/2:(Nz/2-1)); Lz = max(z)-min(z);
[YY,XX] = meshgrid(z,r);
% Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IFFT
ne00 = zeros(Nz, Ns);
ni00 = zeros(Nz, Ns);
phi  = zeros(Nz, Ns);

for it = 1:numel(Ts2D)
    NE_ = Ne00(:,it); NN_ = Ni00(:,it); PH_ = PHI(:,it);
    F_          = (ifft((NE_),Nz));
    ne00(:,it)= real(fftshift(F_));
    F_          = (ifft((NN_),Nz));
    ni00(:,it)  = real(fftshift(F_));
    F_          = (ifft((PH_),Nz));
    phi(:,it) = real(fftshift(F_));
end

% Post processing
ne_00    = zeros(1,Ns);    % Time evolution of ne(r,z) at origin
ni_00    = zeros(1,Ns);    % .
phi_00   = zeros(1,Ns);    % .
E_pot    = zeros(1,Ns);    % Potential energy n^2
E_kin    = zeros(1,Ns);    % Kinetic energy grad(phi)^2
ExB      = zeros(1,Ns);    % ExB drift intensity \propto |\grad \phi|
CFL      = zeros(1,Ns);    % CFL time step
Ddz = 1i*kz;
[~,iz0]  = min(abs(z)); % index of z==0
[~,ikz1] = min(abs(kz-round(1/dkz)*dkz)); % index of kz==1
for it = 1:numel(Ts2D)
    NE_ = squeeze(Ne00(:,it)); NN_ = squeeze(Ni00(:,it)); PH_ = squeeze(PHI(:,it));

    ne_00(it)      = ne00(iz0,it);
    
    ni_00(it)      = ni00(iz0,it);
    
    phi_00(it)     = phi(iz0,it);
    
    E_pot(it)   = pi/Lz*sum(abs(NN_).^2)/Nkz; % integrate through Parseval id
    E_kin(it)   = pi/Lz*sum(abs(Ddz.*PH_).^2)/Nkz;
    ExB(it)     = max(abs(phi(3:end,it)-phi(1:end-2,it))'/(2*dz));
end
E_kin_KZ = abs(Ddz.*PHI(:,it)).^2;
dEdt     = diff(E_pot+E_kin)./diff(Ts2D);

%% Compute growth rate
gamma_Ne = zeros(1,Nkz);
gamma_Ni = zeros(1,Nkz);
gamma_PH = zeros(1,Nkz);
tend   = Ts2D(end); tstart   = 0.6*tend; 
plt = @(x) squeeze(abs(x));
for ikz = 1:Nkz
    gamma_Ne(ikz) = LinearFit_s(Ts2D,plt(Ne00(ikz,:)),tstart,tend);
    gamma_Ni(ikz) = LinearFit_s(Ts2D,plt(Ni00(ikz,:)),tstart,tend);
    gamma_PH(ikz) = LinearFit_s(Ts2D,plt(PHI(ikz,:)),tstart,tend);
end
gamma_Ne = real(gamma_Ne .* (gamma_Ne>=0.0));
gamma_Ni = real(gamma_Ni .* (gamma_Ni>=0.0));
gamma_PH = real(gamma_PH .* (gamma_PH>=0.0));

%% PLOTs2D
if 0
%% Time evolutions
fig = figure; FIGNAME = ['t_evolutions',sprintf('_%.2d',JOBNUM)];
    subplot(221)
        semilogy(Ts2D,abs(ne_00),'-','DisplayName','$n_e^{00}$'); hold on;
        semilogy(Ts2D,abs(ni_00),'-','DisplayName','$n_i^{00}$');
        grid on; xlabel('$t$'); ylabel('$|n_a(x=0,y=0)|$');
    subplot(222)
        semilogy(Ts2D,abs(Ni_gm),'-','DisplayName','$\phi$')
        grid on; xlabel('$t$'); ylabel('$|\tilde n(k_r\approx 1,k_z\approx 1)|$');
    subplot(223)
        semilogy(Ts2D,E_kin+E_pot,'-','DisplayName','$\sum|ik\tilde\phi_i|^2+\sum|\tilde n_i|^2$')
        hold on;
        grid on; xlabel('$t$'); ylabel('$E$'); legend('show');
FMT = '.fig'; save_figure
end
if 1
%% Growth rate
fig = figure; FIGNAME = ['growth_rate',sprintf('_%.2d',JOBNUM)];
plt = @(x) circshift(x,Nkz/2-1);
    subplot(221)
        plot(plt(kz),plt(gamma_Ne),'-'); hold on; xlim([0.0,max(kz)+dkz]);
        grid on; xlabel('$k_z$'); ylabel('$\gamma(N_e^{00})$');
    subplot(222)
        plot(plt(kz),plt(gamma_Ni),'-'); hold on; xlim([0.0,max(kz)+dkz]);
        grid on; xlabel('$k_z$'); ylabel('$\gamma(N_i^{00})$');
    subplot(223)
        plot(plt(kz),plt(gamma_PH),'-'); hold on; xlim([0.0,max(kz)+dkz]);
        grid on; xlabel('$k_z$'); ylabel('$\gamma(\tilde\phi)$'); legend('show');
FMT = '.fig'; save_figure
end
%% GIFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0    = 0;
skip_ = 1; 
DELAY = 0.01*skip_;
FRAMES = floor(t0/dt_samp)+1:skip_:numel(Ts2D);
linestyle = 'o-.';
if 0
%% Density electron
GIFNAME = ['ne',sprintf('_%.2d',JOBNUM)];
FIELD = real(ne00); X = z; T = Ts2D;
FIELDNAME = '$n_e^{00}/\max(n_e^{00})$'; XNAME = '$z$';
XMIN = -L/2-2; XMAX = L/2+1; YMIN = -1.1; YMAX = 1.1;
create_gif_1D
%% Density ion
GIFNAME = ['ni',sprintf('_%.2d',JOBNUM)]; 
FIELD = real(ni00); X = z; T = Ts2D;
FIELDNAME = '$n_i^{00}/\max(n_i^{00})$'; XNAME = '$z$';
XMIN = -L/2-2; XMAX = L/2+1; YMIN = -1.1; YMAX = 1.1;
create_gif_1D
%% Phi
GIFNAME = ['phi',sprintf('_%.2d',JOBNUM)]; 
FIELD = real(phi); X = z; T = Ts2D;
FIELDNAME = '$\phi/\max(\phi)$'; XNAME = '$z$';
XMIN = -L/2-2; XMAX = L/2+1; YMIN = -1.1; YMAX = 1.1;
create_gif_1D
%% Density electron frequency
GIFNAME = ['Ni',sprintf('_%.2d',JOBNUM)]; 
FIELD = (abs(Ni00)); X = (kz); T = Ts2D;
FIELDNAME = '$|N_i^{00}|$'; XNAME = '$k_z$';
XMIN = min(kz)-0.5; XMAX = max(kz)+.5; YMIN = -0.1; YMAX = 1.1;
create_gif_1D
end

if 0
%% Space-Time diagrams at r = 0
plt = @(x) real(x);
fig = figure; FIGNAME = ['z_space_time_diag',sprintf('_%.2d',JOBNUM)];
    [TY,TX] = meshgrid(Ts2D,z);
subplot(221)% density
    pclr = pcolor(TX,TY,(plt(ne00))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$z\,(r=0)$'); ylabel('$t$'); title('$n_e^{00}$');
subplot(222)% density
    pclr = pcolor(TX,TY,(plt(ni00))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$z\,(r=0)$'); ylabel('$t$'); title('$n_i^{00}$');
subplot(223)% density
    pclr = pcolor(TX,TY,(plt(phi))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$z\,(r=0)$'); ylabel('$t$'); title('$\phi$');
FMT = '.fig'; save_figure

%% Space-Time diagrams at kr=0
plt = @(x) abs(x);
fig = figure; FIGNAME = ['kz_space_time_diag',sprintf('_%.2d',JOBNUM)];
    [TY,TX] = meshgrid(Ts2D,kz);
subplot(221)% density
    pclr = pcolor(TX,TY,(plt(Ne00))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$kz$'); ylabel('$t$'); title('$\textrm{Re}(N_e^{00})|_{k_r=0}$');
subplot(222)% density
    pclr = pcolor(TX,TY,(plt(Ni00))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$kz$'); ylabel('$t$'); title('$\textrm{Re}(N_i^{00})|_{k_r=0}$');
subplot(223)% density
    pclr = pcolor(TX,TY,(plt(PH))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$kz$'); ylabel('$t$'); title('$\textrm{Re}(\tilde\phi)|_{k_r=0}$');
FMT = '.fig'; save_figure 
end

if 0
%% Mode time evolution
[~,ik05] = min(abs(kz-0.50));
[~,ik10] = min(abs(kz-1.00));
[~,ik15] = min(abs(kz-1.50));
[~,ik20] = min(abs(kz-2.00));

fig = figure; FIGNAME = ['frame',sprintf('_%.2d',JOBNUM)];
    subplot(221); plt = @(x) abs(squeeze(x));
        semilogy(Ts2D,plt(PHI(ik05,:)),'DisplayName', '$k_z = 0.5$'); hold on
        semilogy(Ts2D,plt(PHI(ik10,:)),'DisplayName', '$k_z = 1.0$')
        semilogy(Ts2D,plt(PHI(ik15,:)),'DisplayName', '$k_z = 1.5$')
        semilogy(Ts2D,plt(PHI(ik20,:)),'DisplayName', '$k_z = 2.0$')
        xlabel('$t$'); ylabel('$|\tilde\phi|$'); legend('show');
        semilogy(Ts2D,plt(PHI(ik05,end)).*exp(gamma_PH(ik05).*(Ts2D-Ts2D(end))),'--k')
        semilogy(Ts2D,plt(PHI(ik10,end)).*exp(gamma_PH(ik10).*(Ts2D-Ts2D(end))),'--k')
        plot(tstart*[1 1],ylim,'k-','LineWidth',0.5);
        plot(tend*[1 1],ylim,'k-','LineWidth',0.5);
    subplot(222); plt = @(x) squeeze(abs(x));
        semilogy(Ts2D,plt(Ni00(ik05,:)),'DisplayName', '$k_z = 0.5$'); hold on
        semilogy(Ts2D,plt(Ni00(ik10,:)),'DisplayName', '$k_z = 1.0$')
        semilogy(Ts2D,plt(Ni00(ik15,:)),'DisplayName', '$k_z = 1.5$')
        semilogy(Ts2D,plt(Ni00(ik20,:)),'DisplayName', '$k_z = 2.0$')        
        xlabel('$t$'); ylabel('$|N_i^{00}|$'); legend('show');
        semilogy(Ts2D,plt(Ni00(ik05,end)).*exp(gamma_Ni(ik05).*(Ts2D-Ts2D(end))),'--k')
        semilogy(Ts2D,plt(Ni00(ik10,end)).*exp(gamma_Ni(ik10).*(Ts2D-Ts2D(end))),'--k')
        plot(tstart*[1 1],ylim,'k-','LineWidth',0.5);
        plot(tend*[1 1],ylim,'k-','LineWidth',0.5);
    subplot(223); plt = @(x) squeeze(abs(x));
        semilogy(Ts2D,plt(Ne00(ik05,:)),'DisplayName', '$k_z = 0.5$'); hold on
        semilogy(Ts2D,plt(Ne00(ik10,:)),'DisplayName', '$k_z = 1.0$')
        semilogy(Ts2D,plt(Ne00(ik15,:)),'DisplayName', '$k_z = 1.5$')
        semilogy(Ts2D,plt(Ne00(ik20,:)),'DisplayName', '$k_z = 2.0$')  
        xlabel('$t$'); ylabel('$|N_e^{00}|$'); legend('show');
        semilogy(Ts2D,plt(Ne00(ik05,end)).*exp(gamma_Ne(ik05).*(Ts2D-Ts2D(end))),'--k')
        semilogy(Ts2D,plt(Ne00(ik10,end)).*exp(gamma_Ne(ik10).*(Ts2D-Ts2D(end))),'--k')
        plot(tstart*[1 1],ylim,'k-','LineWidth',0.5);
        plot(tend*[1 1],ylim,'k-','LineWidth',0.5);
FMT = '.fig'; save_figure
end

if 0
%% Show frame
it = min(70,numel(Ts2D));
fig = figure; FIGNAME = ['frame',sprintf('_%.2d',JOBNUM)];
    subplot(221); plt = @(x) (abs(x));
        plot(kz,plt(PH(:,it)))
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts2D(it))); legend('$\hat\phi$');
    subplot(222); plt = @(x) (abs(x));
        plot(kz,plt(Ni00(:,it)))
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts2D(it))); legend('$\hat n_i^{00}$');
    subplot(223); plt = @(x) (abs(x));
        plot(kz,plt(Ne00(:,it)))
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts2D(it))); legend('$\hat n_e^{00}$');
FMT = '.fig'; save_figure
end