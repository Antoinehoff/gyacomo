addpath(genpath('../matlab')) % ... add
addpath(genpath('../matlab/plots')) % ... add
outfile ='';
%% Directory of the simulation
if 1% Local results
outfile ='';
outfile ='';
% outfile ='artificial_ZF_freeze/sim_A';
% outfile ='simulation_B/cw_FCGK_kp_3.0';
outfile ='nonlin_FCGK/150x75_L_200_P_4_J_2_eta_0.6_nu_1e-01_FCGK_mu_0e+00';
% outfile ='nonlin_PAGK/100x50_L_200_P_4_J_2_eta_0.6_nu_1e-01_PAGK_mu_0e+00';
% outfile ='nonlin_FCGK/100x50_L_200_P_4_J_2_eta_0.6_nu_1e-01_FCGK_mu_0e+00';
% outfile ='simulation_A';
% outfile ='simulation_B/cw_SGGK_like_species';
% outfile ='simulation_A/CO_damping_SGGK';
% outfile ='simulation_A/cw_DGGK_eta_0.5';
% outfile ='simulation_B/cw_DGGK';
    BASIC.RESDIR      = ['../results/',outfile,'/'];
    BASIC.MISCDIR     = ['/misc/HeLaZ_outputs/results/',outfile,'/'];
    system(['mkdir -p ',BASIC.MISCDIR]);
    CMD = ['cp ', BASIC.RESDIR,'outputs* ',BASIC.MISCDIR]; disp(CMD);
    system(CMD);
else% Marconi results
outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/simulation_A/300x150_L_120_P_8_J_4_eta_0.6_nu_1e-01_SGGK_mu_0e+00/out.txt';
% outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/simulation_B/300x150_L_120_P_8_J_4_eta_0.6_nu_5e-01_SGGK_mu_0e+00/out.txt';
BASIC.RESDIR      = ['../',outfile(46:end-8),'/'];
BASIC.MISCDIR     = ['/misc/HeLaZ_outputs/',outfile(46:end-8),'/'];
end

%% Load the results
% Load outputs from jobnummin up to jobnummax
JOBNUMMIN = 00; JOBNUMMAX = 20; 
% JOBNUMMIN = 07; JOBNUMMAX = 20; % For CO damping sim A 
compile_results %Compile the results from first output found to JOBNUMMAX if existing

%% Post-processing
post_processing

%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default_plots_options
disp('Plots')
FMT = '.fig';

if 0
%% Time evolutions and growth rate
plot_time_evolution_and_gr
end

if 1
%% Space time diagramm (fig 11 Ivanov 2020)
TAVG_0 = 1000; TAVG_1 = 2000; % Averaging times duration
plot_radial_transport_and_shear
end

if 0
%% Space time diagramms
cmax = 0.00001 % max of the colorbar for transport
tstart = 0; tend = Ts3D(end); % time window
plot_space_time_diagrams
end

if 0
%% |phi_k|^2 spectra (Kobayashi 2015 fig 3)
% tstart = 0.8*Ts3D(end); tend = Ts3D(end); % Time window
tstart = 2500; tend = 2500;
% tstart = 10000; tend = 12000;
% Chose the field to plot
% FIELD = Ni00;   FNAME = 'Ni00'; FIELDLTX = 'N_i^{00}';
% FIELD = Ne00;   FNAME = 'Ne00'; FIELDLTX = 'N_e^{00}'
FIELD = PHI; FNAME = 'PHI'; FIELDLTX = '\tilde\phi';
% FIELD_ = fft2(Gamma_x); FIELD = FIELD_(1:76,:,:,:); FNAME = 'Gamma_x'; FIELDLTX = '\tilde\Gamma_x';
LOGSCALE = 1; TRENDS = 1; NORMALIZED = 0;
plot_kperp_spectrum
end

if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
t0    =000; iz = 1; ix = 1; iy = 1;
skip_ =4; DELAY = 2e-3*skip_;
[~, it03D] = min(abs(Ts3D-t0)); FRAMES_3D = it03D:skip_:numel(Ts3D);
[~, it05D] = min(abs(Ts5D-t0)); FRAMES_5D = it05D:skip_:numel(Ts5D);
T = Ts3D; FRAMES = FRAMES_3D;
INTERP = 0; NORMALIZED = 0; CONST_CMAP = 0;% Gif options
% Field to plot
% FIELD = dens_e;       NAME = 'ne';   FIELDNAME = 'n_e';
% FIELD = dens_i;       NAME = 'ni';   FIELDNAME = 'n_i';
% FIELD = dens_e-Z_n_e; NAME = 'ne_NZ';FIELDNAME = 'n_e^{NZ}';
FIELD = dens_i-Z_n_i; NAME = 'ni_NZ';FIELDNAME = 'n_i^{NZ}';
% FIELD = temp_e;       NAME = 'Te';   FIELDNAME = 'n_i';
% FIELD = temp_i;       NAME = 'Ti';   FIELDNAME = 'n_i';
% FIELD = temp_e-Z_T_e; NAME = 'Te_NZ';FIELDNAME = 'T_e^{NZ}';
% FIELD = temp_i-Z_T_i; NAME = 'Ti_NZ';FIELDNAME = 'T_i^{NZ}';
% FIELD = ne00;         NAME = 'ne00'; FIELDNAME = 'n_e^{00}';
% FIELD = ni00;   NAME = 'ni00'; FIELDNAME = 'n_i^{00}';
% FIELD = phi;    NAME = 'phi'; FIELDNAME = '\phi';
% FIELD = Z_phi;    NAME = 'Zphi'; FIELDNAME = '\phi_Z';
% FIELD = Gamma_x;  NAME = 'Gamma_x'; FIELDNAME = '\Gamma_x';

% Sliced
% plt = @(x) real(x(ix, :, :,:)); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; 
% plt = @(x) real(x( :,iy, :,:)); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; 
% plt = @(x) real(x( :, :,iz,:)); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; 

% Averaged
% plt = @(x) mean(x,1); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; 
% plt = @(x) mean(x,2); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; 
plt = @(x) mean(x,3); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; 

FIELD = squeeze(plt(FIELD));

% Naming
GIFNAME   = [NAME,sprintf('_%.2d',JOBNUM),'_',PARAMS];

% Create movie (gif or mp4)
create_gif
% create_mov
end

if 0
%% Photomaton : real space

% Chose the field to plot
% FIELD = ni00;          FNAME = 'ni00';    FIELDLTX = 'n_i^{00}';
% FIELD = ne00;          FNAME = 'ne00';    FIELDLTX = 'n_e^{00}'
% FIELD = dens_i;        FNAME = 'ni';      FIELDLTX = 'n_i';
% FIELD = dens_e;        FNAME = 'ne';      FIELDLTX = 'n_e';
% FIELD = dens_e-Z_n_e;   FNAME = 'ne_NZ';   FIELDLTX = 'n_e^{NZ}';
% FIELD = dens_i-Z_n_i;   FNAME = 'ni_NZ';   FIELDLTX = 'n_i^{NZ}';
% FIELD = temp_i;        FNAME = 'Ti';      FIELDLTX = 'T_i';
% FIELD = temp_e;        FNAME = 'Te';      FIELDLTX = 'T_e';
% FIELD = phi;           FNAME = 'phi';     FIELDLTX = '\phi';
% FIELD = Z_phi-phi;     FNAME = 'phi_NZ';  FIELDLTX = '\phi^{NZ}';
FIELD = Z_phi;     FNAME = 'phi_Z';  FIELDLTX = '\phi^{Z}';
% FIELD = Gamma_x;       FNAME = 'Gamma_x'; FIELDLTX = '\Gamma_x';
% FIELD = dens_e-Z_n_e-(Z_phi-phi);       FNAME = 'Non_adiab_part'; FIELDLTX = 'n_e^{NZ}-\phi^{NZ}';

% Chose when to plot it
tf = [9800 10000];

% Sliced
ix = 1; iy = 1; iz = 1;
% plt = @(x,it) real(x(ix, :, :,it)); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; FIELDLTX = [FIELDLTX,'(x=',num2str(round(x(ix))),')']
% plt = @(x,it) real(x( :,iy, :,it)); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; FIELDLTX = [FIELDLTX,'(y=',num2str(round(y(iy))),')']
plt = @(x,it) real(x( :, :,iz,it)); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; FIELDLTX = [FIELDLTX,'(z=',num2str(round(z(iz)/pi)),'\pi)'] 

% Averaged
% plt = @(x,it) mean(x(:,:,:,it),1); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; FIELDLTX = ['\langle ',FIELDLTX,'\rangle_x']
% plt = @(x,it) mean(x(:,:,:,it),2); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; FIELDLTX = ['\langle ',FIELDLTX,'\rangle_y']
% plt = @(x,it) mean(x(:,:,:,it),3); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; FIELDLTX = ['\langle ',FIELDLTX,'\rangle_z'] 


%
TNAME = [];
fig = figure; FIGNAME = [FNAME,TNAME,'_snaps','_',PARAMS]; set(gcf, 'Position',  [100, 100, 1500, 350])
plt_2 = @(x) x;%./max(max(x));
    for i_ = 1:numel(tf)
    [~,it] = min(abs(Ts3D-tf(i_))); TNAME = [TNAME,'_',num2str(Ts3D(it))];
    subplot(1,numel(tf),i_)
        DATA = plt_2(squeeze(plt(FIELD,it)));
        pclr = pcolor((X),(Y),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap(bluewhitered); %caxis([-30,30]);
        xlabel(latexize(XNAME)); ylabel(latexize(YNAME));set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts3D(it)));
    end
    legend(latexize(FIELDLTX));
save_figure
end

if 0
%% Photomaton : k space

% Chose the field to plot
% FIELD = Ni00;   FNAME = 'Ni00'; FIELDLTX = 'N_i^{00}';
% FIELD = Ne00;   FNAME = 'Ne00'; FIELDLTX = 'N_e^{00}'
FIELD = PHI; FNAME = 'PHI'; FIELDLTX = '\tilde\phi';
% FIELD_ = fft2(Gamma_x); FIELD = FIELD_(1:Nx/2+1,:,:,:); FNAME = 'Gamma_x'; FIELDLTX = '\tilde\Gamma_x';
% FIELD_ = fft2(dens_e); FIELD = FIELD_(1:Nx/2+1,:,:,:); FNAME = 'FFT_Dens_e'; FIELDLTX = '\tilde n_e';

% Chose when to plot it
tf = 1000:500:3000;
% tf = 8000;

% Sliced
ix = 1; iy = 1; iz = 1;
% plt = @(x,it) abs(x(ix, :, :,it)); X = KY_YZ; Y = KZ_YZ; XNAME = 'k_y'; YNAME = 'z'; FIELDLTX = [FIELDLTX,'(k_x=',num2str(round(kx(ix))),')'];
% plt = @(x,it) abs(x( :,iy, :,it)); X = KX_XZ; Y = KZ_XZ; XNAME = 'k_x'; YNAME = 'z'; FIELDLTX = [FIELDLTX,'(k_y=',num2str(round(ky(iy))),')'];
plt = @(x,it) abs(x( :, :,iz,it)); X = KX_XY; Y = KY_XY; XNAME = 'k_x'; YNAME = 'k_y'; FIELDLTX = [FIELDLTX,'(z=',num2str((z(iz)/pi)),'\pi)']; 

%
TNAME = [];
fig = figure; FIGNAME = [FNAME,TNAME,'_snaps','_',PARAMS]; set(gcf, 'Position',  [100, 100, 300*numel(tf), 500])
plt_2 = @(x) (fftshift(x,2));
    for i_ = 1:numel(tf)
    [~,it] = min(abs(Ts3D-tf(i_))); TNAME = [TNAME,'_',num2str(Ts3D(it))];
    subplot(1,numel(tf),i_)
        DATA = plt_2(squeeze(plt(FIELD,it)));
        pclr = pcolor(fftshift(X,2),fftshift(Y,2),DATA); set(pclr, 'edgecolor','none');pbaspect([0.5 1 1])
        caxis([0 1]*5e3);
        xlabel(latexize(XNAME)); ylabel(latexize(YNAME));
        if(i_ > 1); set(gca,'ytick',[]); end;
        title(sprintf('$t c_s/R=%.0f$',Ts3D(it)));
    end
    legend(latexize(FIELDLTX));
save_figure
end

if 0
%%
TAVG_0 = 1000; TAVG_1 = 5000; % Averaging times duration
ZF_fourier_analysis
end

if 0
%%
plot_param_evol
end

if 0
%% Radial shear profile
tf = 3000+[900:20:1100];
ymax = 0;
figure
for i_ = 1:numel(tf)
[~,it] = min(abs(Ts3D-tf(i_)));
data = squeeze((mean(dxphi(:,:,1,it),2)));
plot(x,data,'Displayname',sprintf('$t c_s/R=%.0f$',Ts3D(it))); hold on;
ymax = max([ymax abs(min(data)) abs(max(data))]); 
end
xlim([min(x), max(x)]); ylim(1.2*[-1 1]*ymax);
xlabel('$x/\rho_s$'); ylabel('$v_{E\times B,x}$'); grid on
end

if 1
%% zonal vs nonzonal energies for phi(t)
trange = 200:Ns3D;
Ephi_Z           = zeros(1,Ns3D);
Ephi_NZ_kgt0      = zeros(1,Ns3D);    
Ephi_NZ_kgt1      = zeros(1,Ns3D);    
Ephi_NZ_kgt2      = zeros(1,Ns3D);    
high_k_phi       = zeros(1,Ns3D);    
for it = 1:numel(Ts3D)
%     Ephi_NZ(it) = sum(sum(((KY~=0).*abs(PHI(:,:,1,it)).^2)));
%     Ephi_Z(it)  = sum(sum(((KY==0).*abs(PHI(:,:,1,it)).^2)));
    [amp,ikzf] = max(abs((kx~=0).*PHI(:,1,1,it)));
%     Ephi_NZ(it) = sum(sum(((KX~=0).*(KY~=0).*(KX.^2+KY.^2).*abs(PHI(:,:,1,it)).^2)));
    Ephi_NZ_kgt0(it) = sum(sum(((sqrt(KX.^2+KY.^2)>0.0).*(KX~=0).*(KY~=0).*(KX.^2+KY.^2).*abs(PHI(:,:,1,it)).^2)));
    Ephi_NZ_kgt1(it) = sum(sum(((sqrt(KX.^2+KY.^2)>1.0).*(KX~=0).*(KY~=0).*(KX.^2+KY.^2).*abs(PHI(:,:,1,it)).^2)));
    Ephi_NZ_kgt2(it) = sum(sum(((sqrt(KX.^2+KY.^2)>2.0).*(KX~=0).*(KY~=0).*(KX.^2+KY.^2).*abs(PHI(:,:,1,it)).^2)));
%     Ephi_Z(it)  = kx(ikzf)^2*abs(PHI(ikzf,1,1,it)).^2;
    Ephi_Z(it) = sum(sum(((KX~=0).*(KY==0).*(KX.^2).*abs(PHI(:,:,1,it)).^2)));
%     Ephi_NZ(it) = sum(sum(((KX.^2+KY.^2).*abs(PHI(:,:,1,it)).^2)))-Ephi_Z(it);
    high_k_phi(it)  = abs(PHI(18,18,1,it)).^2; % kperp sqrt(2)
%     high_k_phi(it)  = abs(PHI(40,40,1,it)).^2;% kperp 3.5
end
pltx = @(x) x-x(1);
plty = @(x) x./max(squeeze(x));
fig = figure; FIGNAME = ['ZF_turb_energy_vs_time_',PARAMS];
set(gcf, 'Position',  [100, 100, 1400, 500])
subplot(121)
%     yyaxis left
    semilogy(pltx(Ts3D(trange)),plty(Ephi_Z(trange)),'DisplayName',['Zonal, ',CONAME]); hold on;
%     yyaxis right
    semilogy(pltx(Ts3D(trange)),plty(Ephi_NZ_kgt0(trange)),'DisplayName',['NZ, $k_p>0$, ',CONAME]);
    semilogy(pltx(Ts3D(trange)),plty(Ephi_NZ_kgt1(trange)),'DisplayName',['NZ, $k_p>1$, ',CONAME]);
    semilogy(pltx(Ts3D(trange)),plty(Ephi_NZ_kgt2(trange)),'DisplayName',['NZ, $k_p>2$, ',CONAME]);
%     semilogy(pltx(Ts0D),plty(PGAMMA_RI),'DisplayName',['$\Gamma_x$, ',CONAME]);
    title('Energy'); legend('show')
    xlabel('$t c_s/R$'); grid on;% xlim([0 500]);

subplot(122)
    plot(plty(Ephi_Z(trange)),plty(Ephi_NZ_kgt0(trange)));
    title('Phase space'); legend(CONAME)
    xlabel('$E_Z$'); ylabel('$E_{NZ}$'); grid on;% xlim([0 500]);
%     
% subplot(133)
% %     semilogy(pltx(Ts0D),plty(abs(PGAMMA_RI)*SCALE));
%     for ik = [10 20 30]
%     semilogy(pltx(Ts3D(trange)),plty(abs(squeeze(PHI(ik,ik,1,trange))).^2),'DisplayName',[CONAME,', kp=',num2str(sqrt(kx(ik)^2+ky(ik)^2))]); hold on
%     end
%     title('High kperp damping'); legend('show');
%     xlabel('$t c_s/R$'); grid on;% xlim([0 500]);
end

if 0
%% Conservation laws
Nxmax = Nx/2;
Nymax = Ny/2;
mflux_x_i = squeeze(sum((Gamma_x(     1,1:Nxmax,:)+Gamma_x(     1,2:Nxmax+1,:))/2,2)./sum(Gamma_x(     1,1:Nxmax,:)));
mflux_x_o = squeeze(sum((Gamma_x(  Nxmax,1:Nxmax,:)+Gamma_x(  Nxmax,2:Nxmax+1,:))/2,2)./sum(Gamma_x(  Nxmax,1:Nxmax,:)));
mflux_y_i = squeeze(sum((Gamma_y(1:Nxmax,     1,:)+Gamma_y(2:Nxmax+1,     1,:))/2,1)./sum(Gamma_y(1:Nxmax,     1,:)));
mflux_y_o = squeeze(sum((Gamma_y(1:Nxmax,  Nymax,:)+Gamma_y(2:Nxmax+1,  Nymax,:))/2,1)./sum(Gamma_y(1:Nxmax,  Nymax,:)));

mass_cons = mflux_x_i - mflux_x_o + mflux_y_i - mflux_y_o;
%%
figure
plt = @(x) squeeze(mean(mean(x(:,:,1,:),1),2));
subplot(211)
    plot(Ts3D,plt(dens_e+dens_i),'DisplayName','$\delta n_e + \delta n_i$'); hold on;
    plot(Ts3D,plt(ne00+ni00),'DisplayName','$\delta n_e^{00} + \delta n_i^{00}$'); hold on;
    plot(Ts3D,plt(temp_e+temp_i),'DisplayName','$\delta T_e + \delta T_i$'); hold on;
    legend('show'); grid on; xlim([Ts3D(1) Ts3D(end)])
subplot(212);
    plot(Ts3D,mass_cons*(2*pi/Nx/Ny)^2,'DisplayName','in-out'); hold on
%     plot(Ts3D,squeeze(mflux_x_i),'DisplayName','Flux i x');
%     plot(Ts3D,squeeze(mflux_x_o),'DisplayName','Flux o x');
%     plot(Ts3D,squeeze(mflux_y_i),'DisplayName','Flux i y');
%     plot(Ts3D,squeeze(mflux_y_o),'DisplayName','Flux o y');
    legend('show'); grid on; xlim([Ts3D(1) Ts3D(end)]); %ylim([-0.1, 2]*mean(mflux_x_i))
end