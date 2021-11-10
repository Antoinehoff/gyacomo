addpath(genpath('../matlab')) % ... add
addpath(genpath('../matlab/plots')) % ... add
outfile ='';
%% Directory of the simulation
if 1% Local results
outfile ='';
outfile ='Cyclone/100x100x30_5x3_Lx_200_Ly_100_q0_1.4_e_0.18_kN_2.22_kT_6.9_nu_1e-02_DGGK_adiabe';
% outfile ='simulation_A/CO_damping_FCGK';
% outfile ='fluxtube_salphaB_s0/100x100x30_5x3_Lx_200_Ly_100_q0_2.7_e_0.18_kN_2.22_kT_8_nu_1e-01_DGGK_adiabe';
% outfile ='fluxtube_salphaB_s0/50x100x20_5x3_L_300_q0_2.7_e_0.18_kN_2.22_kT_8_nu_1e-01_DGGK_adiabe';
% outfile ='fluxtube_salphaB_s0/50x100x20_5x3_L_300_q0_2.7_e_0.18_kN_2.22_kT_8_nu_1e-01_DGGK_adiabe_Sg';
    LOCALDIR  = ['../results/',outfile,'/'];
    MISCDIR   = ['/misc/HeLaZ_outputs/results/',outfile,'/'];
    system(['mkdir -p ',MISCDIR]);
    CMD = ['rsync ', LOCALDIR,'outputs* ',MISCDIR]; disp(CMD);
    system(CMD);
else% Marconi results
% outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/simulation_A/300x200_L_200_P_8_J_4_eta_0.6_nu_1e-01_PAGK_mu_0e+00/out.txt';
outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/simulation_A/500x500_L_120_P_4_J_2_eta_0.6_nu_1e-01_DGGK_mu_0e+00/out.txt';
% outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/simulation_A/300x150_L_120_P_8_J_4_eta_0.6_nu_1e-01_SGGK_mu_0e+00/out.txt';
% outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/simulation_A/cw_FCGK_kp_3.0/out.txt';
% outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/nonlin_FCGK/150x75_L_200_P_4_J_2_eta_0.6_nu_1e-01_FCGK_mu_0e+00/out.txt';

% BASIC.RESDIR      = ['../',outfile(46:end-8),'/'];
MISCDIR = ['/misc/HeLaZ_outputs/',outfile(46:end-8),'^{NZ}/'];
end

%% Load the results
% Load outputs from jobnummin up to jobnummax
JOBNUMMIN = 00; JOBNUMMAX = 13; 
data = compile_results(MISCDIR,JOBNUMMIN,JOBNUMMAX); %Compile the results from first output found to JOBNUMMAX if existing


%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default_plots_options
disp('Plots')
FMT = '.fig';

if 1
%% Space time diagramm (fig 11 Ivanov 2020)
TAVG_0 = 100; TAVG_1 = 500; % Averaging times duration
fig = plot_radial_transport_and_shear(data,TAVG_0,TAVG_1);
save_figure(data,fig)
end


if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
options.INTERP    = 1;
options.POLARPLOT = 0;
% options.NAME      = '\Gamma_x';
options.NAME      = 'n_i^{NZ}';
options.PLAN      = 'xy';
options.COMP      = 16;
options.TIME      = 0:1:data.Ts3D(end);
% options.TIME      = 140:0.5:160;
data.a = data.EPS * 2000;
create_film(data,options,'.gif')
end

if 0
%% 2D snapshots
% Options
options.INTERP    = 0;
options.POLARPLOT = 0;
% options.NAME      = '\Gamma_x';
options.NAME      = 'n_i^{NZ}';
options.PLAN      = 'kxky';
options.COMP      = 1;
options.TIME      = [50 800 1200];
data.a = data.EPS * 1000;
fig = photomaton(data,options);
save_figure(data,fig)
end

if 0
%% 3D plot on the geometry
options.INTERP    = 1;
options.NAME      = 'n_i';
options.PLANES    = 16;
options.TIME      = 50;
data.rho_o_R      = 1e-3; % Sound larmor radius over Machine size ratio
FIGURE = show_geometry(data,options);
end

if 0
%%
TAVG_0 = 1000; TAVG_1 = 5000; % Averaging times duration
ZF_fourier_analysis
end

if 0
%% Radial shear profile (with moving average)
tf = 1000+[0:100:1000];
ymax = 0;
figure
for i_ = 1:numel(tf)
[~,it] = min(abs(Ts3D-tf(i_)));
data = squeeze((mean(dx2phi(:,:,1,it),2)));
step = 50;
plot(movmean(x,step),movmean(data,step),'Displayname',sprintf('$t c_s/R=%.0f$',Ts3D(it))); hold on;
ymax = max([ymax abs(min(data)) abs(max(data))]); 
end
xlim([min(x), max(x)]); ylim(1.2*[-1 1]*ymax);
xlabel('$x/\rho_s$'); ylabel('$s_{E\times B,x}$'); grid on
end

if 0
%% zonal vs nonzonal energies for phi(t)
t0 = 200;  [~, it0] = min(abs(Ts3D-t0)); 
itend = Ns3D;
trange = it0:itend;
pltx = @(x) x;%-x(1);
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
    title('Energy'); legend('Location','southeast')
    xlim([Ts3D(it0), Ts3D(itend)]);
    ylim([1e-3, 1.5])
    xlabel('$t c_s/R$'); grid on;% xlim([0 500]);

subplot(122)
%     plot(plty(Ephi_Z(trange)),plty(Ephi_NZ_kgt0(trange)));
    plot3(plty(Ephi_Z(trange)),plty(Ephi_NZ_kgt0(trange)),Ts3D(trange));
    title('Phase space'); legend(CONAME)
    xlabel('$E_Z$'); ylabel('$E_{NZ}$'); zlabel('time'); grid on;% xlim([0 500]);
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

if 0
%% Zonal profiles (ky=0)

% Chose the field to plot
FIELD = Ne00.*conj(Ne00);   FNAME = 'Ne002'; FIELDLTX = '|N_e^{00}|^2'
% FIELD = Ni00.*conj(Ni00);   FNAME = 'Ni002'; FIELDLTX = '|N_i^{00}|^2'
% FIELD = abs(PHI); FNAME = 'absPHI'; FIELDLTX = '|\tilde\phi|';
% FIELD = PHI.*conj(PHI); FNAME = 'PHI2'; FIELDLTX = '|\tilde\phi|^2';
% FIELD_ = fft2(Gamma_x); FIELD = FIELD_(1:Nx/2+1,:,:,:); FNAME = 'Gamma_x'; FIELDLTX = '\tilde\Gamma_x';
% FIELD_ = fft2(dens_e); FIELD = FIELD_(1:Nx/2+1,:,:,:); FNAME = 'FFT_Dens_e'; FIELDLTX = '\tilde n_e';

% Chose when to plot it
tf = 200:200:1200;
% tf = 8000;

% Sliced
plt = @(x,it) x( 2:end, 1,1,it)./max(max(x( 2:end, 1,1,:))); XNAME = 'k_x';
%
TNAME = [];
fig = figure; FIGNAME = ['Zonal_',FNAME,TNAME,'_snaps','_',PARAMS]; set(gcf, 'Position',  [100, 100, 600,400])
plt_2 = @(x) (fftshift(x,2));
    for i_ = 1:numel(tf)
    [~,it] = min(abs(Ts3D-tf(i_))); TNAME = [TNAME,'_',num2str(Ts3D(it))];
    DATA = plt_2(squeeze(plt(FIELD,it)));
    semilogy(kx(2:end),DATA,'-','DisplayName',sprintf('$t c_s/R=%.0f$',Ts3D(it))); hold on; grid on;
    xlabel(latexize(XNAME));
    end
title(['$',FIELDLTX,'$ Zonal Spectrum']); legend('show');
save_figure
end

if 0
%% Time evolutions and growth rate
plot_time_evolution_and_gr
end

if 0
%% |phi_k|^2 spectra (Kobayashi 2015 fig 3)
% tstart = 0.8*Ts3D(end); tend = Ts3D(end); % Time window
tstart = 0000; tend = 1000;
% tstart = 10000; tend = 12000;
% Chose the field to plot
% FIELD = Ni00;   FNAME = 'Ni00'; FIELDLTX = 'N_i^{00}';
% FIELD = Ne00;   FNAME = 'Ne00'; FIELDLTX = 'N_e^{00}'
FIELD = PHI; FNAME = 'PHI'; FIELDLTX = '\tilde\phi';
% FIELD_ = fft2(Gamma_x); FIELD = FIELD_(1:76,:,:,:); FNAME = 'Gamma_x'; FIELDLTX = '\tilde\Gamma_x';
LOGSCALE = 1; TRENDS = 1; NORMALIZED = 0;
plot_kperp_spectrum
end