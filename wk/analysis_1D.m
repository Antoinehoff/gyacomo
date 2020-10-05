default_plots_options
%% load results
JOBNUM = 00;
filename = [BASIC.SIMID,'_','%.2d.h5'];
filename = sprintf(filename,JOBNUM); disp(['Analysing ',filename])
[Nipj, p_, j_, kr, kz, Ts, dt] = load_5D_data(filename, 'moments_i');
Nepj                           = load_5D_data(filename, 'moments_e');
Ni      = squeeze(Nipj(1,1,:,:,:));
Ne      = squeeze(Nepj(1,1,:,:,:));
PH      = squeeze(load_2D_data(filename, 'phi'));
Ts      = Ts';
Ns      = numel(Ts);
dt_samp = mean(diff(Ts));
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
ne   = zeros(Nz, Ns);
ni   = zeros(Nz, Ns);
phi  = zeros(Nz, Ns);

for it = 1:numel(PH(1,:))
    NE_ = Ne(:,it); NN_ = Ni(:,it); PH_ = PH(:,it);
    F_          = (ifft((NE_),Nz));
    ne(:,it)= real(fftshift(F_));
    F_          = (ifft((NN_),Nz));
    ni(:,it)  = real(fftshift(F_));
    F_          = (ifft((PH_),Nz));
    phi(:,it) = real(fftshift(F_));
end

% Post processing
ne_00    = zeros(1,Ns);    % Time evolution of ne(r,z) at origin
Ne_gm    = zeros(1,Ns);    % Time evolution of Ne(k) max gamma (max real)
ni_00    = zeros(1,Ns);    % .
Ni_gm    = zeros(1,Ns);    % .
phi_00   = zeros(1,Ns);    % .
E_pot    = zeros(1,Ns);    % Potential energy n^2
E_kin    = zeros(1,Ns);    % Kinetic energy grad(phi)^2
ExB      = zeros(1,Ns);    % ExB drift intensity \propto |\grad \phi|
CFL      = zeros(1,Ns);    % CFL time step
Ddz = 1i*kz;
[~,iz0]  = min(abs(z)); % index of z==0
[~,ikz1] = min(abs(kz-round(1/dkz)*dkz)); % index of kz==1
for it = 1:numel(PH(1,:))
    NE_ = squeeze(Ne(:,it)); NN_ = squeeze(Ni(:,it)); PH_ = squeeze(PH(:,it));

    ne_00(it)      = ne(iz0,it);
    Ne_gm(it)      = max(real(Ne(:,it)));
    
    ni_00(it)      = ni(iz0,it);
    Ni_gm(it)      = max(real(Ni(:,it)));
    
    phi_00(it)     = phi(iz0,it);
    
    E_pot(it)   = pi/Lz*sum(abs(NN_).^2)/Nkz; % integrate through Parseval id
    E_kin(it)   = pi/Lz*sum(abs(Ddz.*PH_).^2)/Nkz;
    ExB(it)     = max(abs(phi(3:end,it)-phi(1:end-2,it))'/(2*dz));
end
E_kin_KZ = abs(Ddz.*PH(:,it)).^2;
dEdt     = diff(E_pot+E_kin)./diff(Ts);

% Growth rate
gamma_Ne = zeros(1,Nkz);
gamma_Ni = zeros(1,Nkz);
gamma_PH = zeros(1,Nkz);
for ikz = 1:Nkz
    gamma_Ne(ikz) = real(LinearFit_s(Ts,Ne(ikz,:)));
    gamma_Ni(ikz) = real(LinearFit_s(Ts,Ni(ikz,:)));
    gamma_PH(ikz) = real(LinearFit_s(Ts,PH(ikz,:)));
end

%% PLOTS
%% Time evolutions
fig = figure; FIGNAME = ['t_evolutions',sprintf('_%.2d',JOBNUM)];
    subplot(221)
        semilogy(Ts,abs(ne_00),'-','DisplayName','$n_e^{00}$'); hold on;
        semilogy(Ts,abs(ni_00),'-','DisplayName','$n_i^{00}$');
        grid on; xlabel('$t$'); ylabel('$|n_a(x=0,y=0)|$');
    subplot(222)
        semilogy(Ts,abs(Ni_gm),'-','DisplayName','$\phi$')
        grid on; xlabel('$t$'); ylabel('$|\tilde n(k_r\approx 1,k_z\approx 1)|$');
    subplot(223)
        semilogy(Ts,E_kin+E_pot,'-','DisplayName','$\sum|ik\tilde\phi_i|^2+\sum|\tilde n_i|^2$')
        hold on;
        grid on; xlabel('$t$'); ylabel('$E$'); legend('show');
FMT = '.fig'; save_figure

%% Growth rate
fig = figure; FIGNAME = ['growth_rate',sprintf('_%.2d',JOBNUM)];
    subplot(221)
        plot(kz,gamma_Ne,'-'); hold on;
        grid on; xlabel('$k_z$'); ylabel('$\gamma(N_e^{00})$');
    subplot(222)
        plot(kz,gamma_Ni,'-'); hold on;
        grid on; xlabel('$k_z$'); ylabel('$\gamma(N_i^{00})$');
    subplot(223)
        plot(kz,gamma_PH,'-'); hold on;
        grid on; xlabel('$k_z$'); ylabel('$\gamma(\tilde\phi)$'); legend('show');
FMT = '.fig'; save_figure

%% GIFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DELAY = 0.07; skip_ = 10;
FRAMES = 1:skip_:numel(Ts);
if 0
%% Density electron
GIFNAME = ['ne',sprintf('_%.2d',JOBNUM)];
FIELD = real(ne); X = z; T = Ts;
FIELDNAME = '$n_e^{00}$'; XNAME = '$z$';
XMIN = -20; XMAX = 20; YMIN = -0.6; YMAX = 1.1;
create_gif_1D
%% Density ion
GIFNAME = ['ni',sprintf('_%.2d',JOBNUM)]; 
FIELD = real(ni); X = z; T = Ts;
FIELDNAME = '$n_i^{00}$'; XNAME = '$z$';
XMIN = -20; XMAX = 20; YMIN = -0.6; YMAX = 1.1;
create_gif_1D
%% Phi
GIFNAME = ['phi',sprintf('_%.2d',JOBNUM)]; 
FIELD = real(phi); X = z; T = Ts;
FIELDNAME = '$\phi$'; XNAME = '$z$';
XMIN = -20; XMAX = 20; YMIN = -1.1; YMAX = 1.1;
create_gif_1D
%% Density electron frequency
GIFNAME = ['Ni',sprintf('_%.2d',JOBNUM)]; 
FIELD = fftshift(real(Ni)); X = kr\z; T = Ts;
FIELDNAME = '$N_i^{00}$'; XNAME = '$k_z$';
create_gif_1D
end

if 0
%% Space-Time diagrams at r = 0
plt = @(x) real(x);
fig = figure; FIGNAME = ['z_space_time_diag',sprintf('_%.2d',JOBNUM)];
    [TY,TX] = meshgrid(Ts,z);
subplot(221)% density
    pclr = pcolor(TX,TY,(plt(ne))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$z\,(r=0)$'); ylabel('$t$'); title('$n_e^{00}$');
subplot(222)% density
    pclr = pcolor(TX,TY,(plt(ni))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$z\,(r=0)$'); ylabel('$t$'); title('$n_i^{00}$');
subplot(223)% density
    pclr = pcolor(TX,TY,(plt(phi))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$z\,(r=0)$'); ylabel('$t$'); title('$\phi$');
FMT = '.fig'; save_figure

%% Space-Time diagrams at kr=0
plt = @(x) abs(x);
fig = figure; FIGNAME = ['kz_space_time_diag',sprintf('_%.2d',JOBNUM)];
    [TY,TX] = meshgrid(Ts,kz);
subplot(221)% density
    pclr = pcolor(TX,TY,(plt(Ne))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$kz$'); ylabel('$t$'); title('$\textrm{Re}(N_e^{00})|_{k_r=0}$');
subplot(222)% density
    pclr = pcolor(TX,TY,(plt(Ni))); set(pclr, 'edgecolor','none'); colorbar;
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
    subplot(221); plt = @(x) (abs(x));
        semilogy(Ts,plt(PH(ik05,:)),'DisplayName', '$k_z = 0.5$'); hold on
        semilogy(Ts,plt(PH(ik10,:)),'DisplayName', '$k_z = 1.0$')
        semilogy(Ts,plt(PH(ik15,:)),'DisplayName', '$k_z = 1.5$')
        semilogy(Ts,plt(PH(ik20,:)),'DisplayName', '$k_z = 2.0$')
        xlabel('$t$'); ylabel('$\hat\phi$'); legend('show');
        semilogy(Ts,abs(PH(ik05,end)).*exp(gamma_PH(ik05).*(Ts-Ts(end))),'--k')
        semilogy(Ts,abs(PH(ik10,end)).*exp(gamma_PH(ik10).*(Ts-Ts(end))),'--k')
    subplot(222); plt = @(x) (abs(x));
        semilogy(Ts,plt(Ni(ik05,:)),'DisplayName', '$k_z = 0.5$'); hold on
        semilogy(Ts,plt(Ni(ik10,:)),'DisplayName', '$k_z = 1.0$')
        semilogy(Ts,plt(Ni(ik15,:)),'DisplayName', '$k_z = 1.5$')
        semilogy(Ts,plt(Ni(ik20,:)),'DisplayName', '$k_z = 2.0$')        
        xlabel('$t$'); ylabel('$\hat n_i^{00}$'); legend('show');
    subplot(223); plt = @(x) (abs(x));
        semilogy(Ts,plt(Ne(ik05,:)),'DisplayName', '$k_z = 0.5$'); hold on
        semilogy(Ts,plt(Ne(ik10,:)),'DisplayName', '$k_z = 1.0$')
        semilogy(Ts,plt(Ne(ik15,:)),'DisplayName', '$k_z = 1.5$')
        semilogy(Ts,plt(Ne(ik20,:)),'DisplayName', '$k_z = 2.0$')  
        xlabel('$t$'); ylabel('$\hat n_e^{00}$'); legend('show');
FMT = '.fig'; save_figure
end

if 0
%% Show frame
it = min(70,numel(Ts));
fig = figure; FIGNAME = ['frame',sprintf('_%.2d',JOBNUM)];
    subplot(221); plt = @(x) (abs(x));
        plot(kz,plt(PH(:,it)))
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts(it))); legend('$\hat\phi$');
    subplot(222); plt = @(x) (abs(x));
        plot(kz,plt(Ni(:,it)))
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts(it))); legend('$\hat n_i^{00}$');
    subplot(223); plt = @(x) (abs(x));
        plot(kz,plt(Ne(:,it)))
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts(it))); legend('$\hat n_e^{00}$');
FMT = '.fig'; save_figure
end