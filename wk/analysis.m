%% load results
JOBNUM = 00;
filename = [BASIC.SIMID,'_','%.2d.h5'];
filename = sprintf(filename,JOBNUM); disp(['Analysing ',filename])
[Nipj, p_, j_, kr, kz, Ts, dt] = load_5D_data(filename, 'moments_i');
Nepj                           = load_5D_data(filename, 'moments_e');
Ni      = squeeze(Nipj(1,1,:,:,:));
Ne      = squeeze(Nepj(1,1,:,:,:));
PH      = load_2D_data(filename, 'phi');
Ts      = Ts';
Ns      = numel(Ts);
dt_samp = mean(diff(Ts));
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
Ne_ST_kr = zeros(Nkr,Ns);   % Space-Time of max_kz Ne(k)
Ne_ST_kz = zeros(Nkz,Ns);   % ''            max_kr Ne(k)
ne_ST_r  = zeros(Nr,Ns);   % Space-Time of ne(z==0)
ne_ST_z  = zeros(Nz,Ns);   % ''               r==0
ne_00    = zeros(1,Ns);    % Time evolution of ne(r,z) at origin
Ne_gm    = zeros(1,Ns);    % Time evolution of Ne(k) max gamma (max real)
Ni_ST_kr = zeros(Nkr,Ns);   % same for ions
Ni_ST_kz = zeros(Nkz,Ns);   % .
ni_ST_r  = zeros(Nr,Ns);   % . 
ni_ST_z  = zeros(Nz,Ns);   % .
ni_00    = zeros(1,Ns);    % .
Ni_gm    = zeros(1,Ns);    % .
PH_ST_kr = zeros(Nkr,Ns);   % same for ES-potential
PH_ST_kz = zeros(Nkz,Ns);   % . 
phi_ST_r = zeros(Nr,Ns);   % .
phi_ST_z = zeros(Nz,Ns);   % .
phi_00   = zeros(1,Ns);    % .
Sne_norm = zeros(1,Ns);    % Time evolution of the amp of e nonlin term
Sni_norm = zeros(1,Ns);    % Time evolution of the amp of i nonlin term
E_pot    = zeros(1,Ns);    % Potential energy n^2
E_kin    = zeros(1,Ns);    % Kinetic energy grad(phi)^2
ExB      = zeros(1,Ns);    % ExB drift intensity \propto |\grad \phi|
CFL      = zeros(1,Ns);    % CFL time step
Ddr = 1i*KR; Ddz = 1i*KZ; lapl   = Ddr.^2 + Ddz.^2; 
[~,ir0]  = min(abs(r)); % index of r==0
[~,iz0]  = min(abs(z)); % index of z==0
[~,ikr1] = min(abs(kr-round(1/dkr)*dkr)); % index of kr==1 
[~,ikz1] = min(abs(kz-round(1/dkz)*dkz)); % index of kz==1
for it = 1:numel(PH(1,1,:))
    NE_ = Ne(:,:,it); NN_ = Ni(:,:,it); PH_ = PH(:,:,it);

    ne_ST_r  (:,it)= ne(:,iz0,it);  ne_ST_z  (:,it)= ne(ir0,:,it);
    Ne_ST_kr (:,it)= max(real(Ne(:,:,it)),[],2); 
    Ne_ST_kz (:,it)= max(real(Ne(:,:,it)),[],1);
    ne_00(it)      = ne(ir0,iz0,it);
    Ne_gm(it)      = max(max(real(Ne(:,:,it))));
    
    ni_ST_r  (:,it)= ni(:,iz0,it);  ni_ST_z  (:,it)= ni(ir0,:,it);
    Ni_ST_kr (:,it)= max(real(Ni(:,:,it)),[],2); 
    Ni_ST_kz (:,it)= max(real(Ni(:,:,it)),[],1);
    ni_00(it)      = ni(ir0,iz0,it);
    Ni_gm(it)      = max(max(real(Ni(:,:,it))));
    
    phi_ST_r(:,it) = phi(:,iz0,it); phi_ST_z(:,it) = phi(ir0,:,it);
    PH_ST_kr(:,it) = max(real(PH(:,:,it)),[],2); 
    PH_ST_kz(:,it) = max(real(PH(:,:,it)),[],1);
    phi_00(it)     = phi(ir0,iz0,it);
    
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
        semilogy(Ts,abs(Ni_gm),'-','DisplayName','$\phi$')
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
% fig = figure; FIGNAME = ['Energy_kin_KZ',sprintf('_%.2d',JOBNUM)];
%     semilogy(kr(floor(end/2)+1:end),E_kin_KR(floor(end/2)+1:end),'o','DisplayName','$\sum_y\langle|ik\tilde\phi_i|^2\rangle_t$')
%     hold on;
%     loglog(kz(floor(end/2)+1:end),E_kin_KZ(floor(end/2)+1:end),'o','DisplayName','$\sum_x\langle|ik\tilde\phi_i|^2\rangle_t$')
%     grid on; xlabel('$k$');  legend('show');
% FMT = '.fig'; save_figure

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

%% Space-Time diagrams for max_kz(Real)
plt = @(x) (x);
fig = figure; FIGNAME = ['kr_space_time_diag',sprintf('_%.2d',JOBNUM)];
    [TY,TX] = meshgrid(Ts,kr);
subplot(221)% density
    pclr = pcolor(TX,TY,(plt(Ne_ST_kr))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$kr$'); ylabel('$t$'); title('$\max_{k_z}(\textrm{Re}(N_e^{00}))$');
subplot(222)% density
    pclr = pcolor(TX,TY,(plt(Ni_ST_kr))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$kr$'); ylabel('$t$'); title('$\max_{k_z}(\textrm{Re}(N_i^{00}))$');
subplot(223)% density
    pclr = pcolor(TX,TY,(plt(PH_ST_kr))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$kr$'); ylabel('$t$'); title('$\max_{k_z}(\textrm{Re}(\tilde\phi$))');
FMT = '.fig'; save_figure

%% Space-Time diagrams at max_kr(Real)
plt = @(x) (x);
fig = figure; FIGNAME = ['kz_space_time_diag',sprintf('_%.2d',JOBNUM)];
    [TY,TX] = meshgrid(Ts,kz);
subplot(221)% density
    pclr = pcolor(TX,TY,(plt(Ne_ST_kz))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$kz$'); ylabel('$t$'); title('$\max_{k_r}(\textrm{Re}(N_e^{00}))$');
subplot(222)% density
    pclr = pcolor(TX,TY,(plt(Ni_ST_kz))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$kz$'); ylabel('$t$'); title('$\max_{k_r}(\textrm{Re}(N_i^{00}))$');
subplot(223)% density
    pclr = pcolor(TX,TY,(plt(PH_ST_kz))); set(pclr, 'edgecolor','none'); colorbar;
    xlabel('$kz$'); ylabel('$t$'); title('$\max_{k_r}(\textrm{Re}(\tilde\phi))$');
FMT = '.fig'; save_figure


if 0
%% Show frame
it = min(70,numel(Ts));
fig = figure; FIGNAME = ['frame',sprintf('_%.2d',JOBNUM)];
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
DELAY = 0.07; skip_ = 10;
FRAMES = 1:skip_:numel(Ts);
if 0
%% GIFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Density electron
GIFNAME = ['ne',sprintf('_%.2d',JOBNUM)];
FIELD = real(ne); X = XX; Y = YY; T = Ts;
FIELDNAME = '$n_e^{00}$'; XNAME = '$r$'; YNAME = '$z$';
create_gif
%% Density ion
GIFNAME = ['ni',sprintf('_%.2d',JOBNUM)]; 
FIELD = real(ni); X = XX; Y = YY; T = Ts;
FIELDNAME = '$n_i^{00}$'; XNAME = '$r$'; YNAME = '$z$';
create_gif
%% Phi
GIFNAME = ['phi',sprintf('_%.2d',JOBNUM)]; 
FIELD = real(phi); X = XX; Y = YY; T = Ts;
FIELDNAME = '$\phi$'; XNAME = '$r$'; YNAME = '$z$';
create_gif
%% Density electron frequency
GIFNAME = ['Ni',sprintf('_%.2d',JOBNUM)]; 
FIELD = fftshift(real(Ni)); X = fftshift(KR); Y = fftshift(KZ); T = Ts;
FIELDNAME = '$N_i^{00}$'; XNAME = '$k_r$'; YNAME = '$k_z$';
create_gif
end