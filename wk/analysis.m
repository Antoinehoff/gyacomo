%% load results
JOBNUM = 00;
filename = [BASIC.SIMID,'_','%.2d.h5'];
filename = sprintf(filename,JOBNUM); disp(['Analysing ',filename])
[Nipj, p_, j_, kr, kz, Ts] = load_5D_data(filename, 'moments_i');
Nepj                       = load_5D_data(filename, 'moments_e');
Ni      = squeeze(Nipj(1,1,:,:,:));
Ne      = squeeze(Nepj(1,1,:,:,:));
PH      = load_2D_data(filename, 'phi');
Ts      = Ts';
Ns      = numel(Ts);
dt      = mean(diff(Ts));
if strcmp(OUTPUTS.write_non_lin,'.true.')
    Sipj    = load_5D_data(filename, 'Sipj');
    Sepj    = load_5D_data(filename, 'Sepj');
    SNi      = squeeze(Sipj(1,1,:,:,:));
    SNe      = squeeze(Sepj(1,1,:,:,:));
end
%% Build grids
Nkr = numel(kr); Nkz = numel(kz);
[KZ,KR] = meshgrid(kz,kr);
Lkr = max(kr)-min(kr); Lkz = max(kz)-min(kz);
dkr = Lkr/(Nkr-1); dkz = Lkz/(Nkz-1);
Lk = max(Lkr,Lkz);

dr = 2*pi/Lk; dz = 2*pi/Lk;
Nr = max(Nkr,Nkz);         Nz = Nr;
r = dr*(-Nr/2:(Nr/2-1)); Lr = max(r)-min(r);
z = dz*(-Nz/2:(Nz/2-1)); Lz = max(z)-min(z);
[YY,XX] = meshgrid(z,r);
%% Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IFFT
ne   = zeros(Nr,Nz);
ni     = zeros(Nr,Nz);
phi    = zeros(Nr,Nz);

for it = 1:numel(PH(1,1,:))
    NE_ = Ne(:,:,it); NN_ = Ni(:,:,it); PH_ = PH(:,:,it);
    F_          = (ifft2((NE_),Nr,Nz));
    ne(:,:,it)= real(fftshift(F_));
    F_          = (ifft2((NN_),Nr,Nz));
    ni(:,:,it)  = real(fftshift(F_));
    F_          = (ifft2((PH_),Nr,Nz));
    phi(:,:,it) = real(fftshift(F_));
end

%% Post processing
phi_ST_r = zeros(Nr,Ns);   % Space-Time diagram of ES potential
ne_ST_r  = zeros(Nr,Ns);   % Space-Time diagram of density
ni_ST_r  = zeros(Nr,Ns);   % Space-Time diagram of density
phi_ST_z = zeros(Nz,Ns);   % Space-Time diagram of ES potential
ne_ST_z  = zeros(Nz,Ns);   % Space-Time diagram of density
ni_ST_z  = zeros(Nz,Ns);   % Space-Time diagram of density
phi_00 = zeros(1,Ns);    % Time evolution of ES potential at origin
ne_00  = zeros(1,Ns);    % Time evolution of density at origin
ni_00  = zeros(1,Ns);    % Time evolution of density at origin
[~,ir0] = min(abs(r)); [~,iz0] = min(abs(z));
Ne_11  = zeros(1,Ns);    % Time evolution of F density at 1,1
Ni_11  = zeros(1,Ns);    % Time evolution of F density at 1,1
[~,ikr1] = min(abs(kr-round(1/dkr)*dkr)); [~,ikz1] = min(abs(kz-round(1/dkz)*dkz));
Sni_norm= zeros(1,Ns);   % Time evolution of the amp of density nonlin term
Sne_norm= zeros(1,Ns);   % Time evolution of the amp of vorti. nonlin term
E_pot  = zeros(1,Ns);    % Potential energy n^2
E_kin  = zeros(1,Ns);    % Kinetic energy grad(phi)^2
ExB    = zeros(1,Ns);    % ExB drift intensity \propto |\grad \phi|
CFL    = zeros(1,Ns);    % CFL time step
Ddr = 1i*KR; Ddz = 1i*KZ; lapl   = Ddr.^2 + Ddz.^2; 
for it = 1:numel(PH(1,1,:))
    NE_ = Ne(:,:,it); NN_ = Ni(:,:,it); PH_ = PH(:,:,it);
    phi_ST_r(:,it) = phi(:,iz0,it); phi_ST_z(:,it) = phi(ir0,:,it);
    ne_ST_r  (:,it)= ne(:,iz0,it);  ne_ST_z  (:,it)= ne(ir0,:,it);
    ni_ST_r  (:,it)= ni(:,iz0,it);  ni_ST_z  (:,it)= ni(ir0,:,it);
    phi_00(it)   = phi(ir0,iz0,it);
    ne_00(it)    = ne(ir0,iz0,it);
    ni_00(it)    = ni(ir0,iz0,it);
    Ne_11(it)    = Ne(ikr1,ikz1,it);
    Ni_11(it)    = Ni(ikr1,ikz1,it);
if strcmp(OUTPUTS.write_non_lin,'.true.')
    Sni_norm(it) = sum(sum(abs(SNi(:,:,it))));
    Sne_norm(it) = sum(sum(abs(SNe(:,:,it))));
end
    E_pot(it)   = pi/Lr/Lz*sum(sum(abs(NN_).^2))/Nkr/Nkr; % integrate through Parseval id
    E_kin(it)   = pi/Lr/Lz*sum(sum(abs(Ddr.*PH_).^2+abs(Ddz.*PH_).^2))/Nkr/Nkr;
    ExB(it)     = max(max(max(abs(phi(3:end,:,it)-phi(1:end-2,:,it))/(2*dr))),max(max(abs(phi(:,3:end,it)-phi(:,1:end-2,it))'/(2*dz))));
end
E_kin_KZ = mean(mean(abs(Ddr.*PH(:,:,it)).^2+abs(Ddz.*PH(:,:,it)).^2,3),1);
E_kin_KR = mean(mean(abs(Ddr.*PH(:,:,it)).^2+abs(Ddz.*PH(:,:,it)).^2,3),2);
dEdt     = diff(E_pot+E_kin)./diff(Ts);
%% PLOTS
%% Time evolutions
fig = figure; FIGNAME = ['t_evolutions',sprintf('_%.2d',JOBNUM)];
    subplot(221)
        semilogy(Ts,abs(ne_00),'-','DisplayName','$n_e^{00}$'); hold on;
        semilogy(Ts,abs(ni_00),'-','DisplayName','$n_i^{00}$');
        grid on; xlabel('$t$'); ylabel('$|n_a(x=0,y=0)|$');
    subplot(222)
        semilogy(Ts,abs(Ni_11),'-','DisplayName','$\phi$')
        grid on; xlabel('$t$'); ylabel('$|\tilde n(k_r\approx 1,k_z\approx 1)|$');
    subplot(223)
        semilogy(Ts,E_kin+E_pot,'-','DisplayName','$\sum|ik\tilde\phi_i|^2+\sum|\tilde n_i|^2$')
        hold on;
        grid on; xlabel('$t$'); ylabel('$E$'); legend('show');
if strcmp(OUTPUTS.write_non_lin,'.true.')
    subplot(224)
        semilogy(Ts,Sne_norm,'-','DisplayName','$\sum|S_e^{00}|$'); 
        hold on;
        semilogy(Ts,Sni_norm,'-','DisplayName','$\sum|S_i^{00}|$');
        grid on; xlabel('$t$'); ylabel('$S$'); legend('show');
end
FMT = '.fig'; save_figure

%% Spectra energy
fig = figure; FIGNAME = ['Energy_kin_KZ',sprintf('_%.2d',JOBNUM)];
    semilogy(kr(floor(end/2)+1:end),E_kin_KR(floor(end/2)+1:end),'o','DisplayName','$\sum_y\langle|ik\tilde\phi_i|^2\rangle_t$')
    hold on;
    loglog(kz(floor(end/2)+1:end),E_kin_KZ(floor(end/2)+1:end),'o','DisplayName','$\sum_x\langle|ik\tilde\phi_i|^2\rangle_t$')
    grid on; xlabel('$k$');  legend('show');
FMT = '.fig'; save_figure

%% CFL condition
fig = figure; FIGNAME = ['CFL',sprintf('_%.2d',JOBNUM)];
    semilogy(Ts,dz./ExB,'-','DisplayName','$|\nabla \phi|\Delta y$');
    hold on;
    plot(Ts,dt*ones(1,numel(Ts)),'--k','DisplayName','$\Delta t$');
    grid on; xlabel('$t$'); ylabel('$\Delta t$'); legend('show');
FMT = '.fig'; save_figure

%% Space-Time diagrams at z = 0
plt = @(x) real(x);
fig = figure; FIGNAME = ['r_space_time_diag',sprintf('_%.2d',JOBNUM)];
    [TY,TX] = meshgrid(Ts,z);
subplot(221)% density
    pclr = pcolor(TX,TY,(plt(ne_ST_r))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$r\,(z=0)$'); ylabel('$t$'); title('$n_e^{00}$');
subplot(222)% density
    pclr = pcolor(TX,TY,(plt(ni_ST_r))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$r\,(z=0)$'); ylabel('$t$'); title('$n_i^{00}$');
subplot(223)% density
    pclr = pcolor(TX,TY,(plt(phi_ST_r))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$r\,(z=0)$'); ylabel('$t$'); title('$\phi$');
FMT = '.fig'; save_figure

%% Space-Time diagrams at r = 0
plt = @(x) real(x);
fig = figure; FIGNAME = ['z_space_time_diag',sprintf('_%.2d',JOBNUM)];
    [TY,TX] = meshgrid(Ts,r);
subplot(221)% density
    pclr = pcolor(TX,TY,(plt(ne_ST_z))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$z\,(r=0)$'); ylabel('$t$'); title('$n_e^{00}$');
subplot(222)% density
    pclr = pcolor(TX,TY,(plt(ni_ST_z))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$z\,(r=0)$'); ylabel('$t$'); title('$n_i^{00}$');
subplot(223)% density
    pclr = pcolor(TX,TY,(plt(phi_ST_z))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$z\,(r=0)$'); ylabel('$t$'); title('$\phi$');
FMT = '.fig'; save_figure

%% phi
fig = figure; FIGNAME = ['phi_ST',sprintf('_%.2d',JOBNUM)];
    [TY,TX] = meshgrid(Ts,z);
    pclr = pcolor(TX,TY,(plt(phi_ST_r))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$x\,(y=0)$'); ylabel('$t$'); title('$\phi$');
FMT = '.fig'; save_figure

if 0
%% Show frame
it = min(1,numel(Ts));
fig = figure; FIGNAME = ['frame',sprintf('_%.2d',SID)];
    subplot(221); plt = @(x) fftshift((real(x)));
        pclr = pcolor(fftshift(KR),fftshift(KZ),plt(PH(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts(it))); legend('$\hat\phi$');
    subplot(222); plt = @(x) fftshift(real(x));
        pclr = pcolor(fftshift(KR),fftshift(KZ),plt(Ni(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts(it))); legend('$\hat n_i^{00}$');
    subplot(223); plt = @(x) fftshift(real(x));
        pclr = pcolor(fftshift(KR),fftshift(KZ),plt(Ne(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('t=%.3d',Ts(it))); legend('$\hat n_e^{00}$');
if strcmp(OUTPUTS.write_non_lin,'.true.')
    subplot(224); plt = @(x) fftshift(log10(abs(x)));
        pclr = pcolor(fftshift(KR),fftshift(KZ),plt(SNi(:,:,it))); set(pclr, 'edgecolor','none'); colorbar;
        xlabel('$k_r$'); ylabel('$k_z$');legend('$\hat S_i^{00}$');
end
FMT = '.fig'; save_figure
end
%%
DELAY = 0.1; skip_ = 2;
if 0
%% GIFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Density electron
GIFNAME = ['ne',sprintf('_%.2d',JOBNUM)]; FIELDNAME = '$n_e^{00}$';
FIELD = real(ne(:,:,1:skip_:end)); X = XX; Y = YY; T = Ts(1:skip_:end);
create_gif
%% Density ion
GIFNAME = ['ni',sprintf('_%.2d',JOBNUM)]; FIELDNAME = '$n_i^{00}$';
FIELD = real(ni(:,:,1:skip_:end)); X = XX; Y = YY; T = Ts(1:skip_:end);
create_gif
%% Phi
GIFNAME = ['phi',sprintf('_%.2d',JOBNUM)]; FIELDNAME = '$\phi$';
FIELD = real(phi(:,:,1:skip_:end)); X = XX; Y = YY; T = Ts(1:skip_:end);
create_gif
end