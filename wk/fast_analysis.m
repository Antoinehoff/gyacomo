gyacomodir = pwd; gyacomodir = gyacomodir(1:end-2); % get code directory
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add
default_plots_options
% Partition of the computer where the data have to be searched
% PARTITION='/Users/ahoffmann/gyacomo/results/';
% PARTITION='/home/ahoffman/gyacomo/results/';
% PARTITION='/misc/gyacomo23_outputs/';
% resdir = 'AE_3x2x128x32x24/PT';
% resdir = 'AE_3x2x128x32x24/PT';
% resdir = 'AE_5x3x128x32x24/NT';
% resdir = 'IS_5x3x128x32x24/NT';
% PARTITION='/home/ahoffman/gyacomo/simulations/';
% resdir ='cheap_CBC_baseline';
% resdir ='test_HEL_closure';
% resdir ='dmax_closure/';
% PARTITION = '/misc/gyacomo23_outputs/reduced_fluid_paper/';
% resdir = '/Npol_study/AE_CBC_s0_beta0_P4J2';
% resdir = '/Npol_study/RF_CBC_s0_beta0/Npol_11';
% Triangularity paper

PARTITION = '/misc/gyacomo23_outputs/triangularity_paper/';
% Nominal parameters
% resdir = 'ion_scale/3x2x256x64x32/0T';
% resdir = 'ion_scale/5x3x256x64x32/0T';
resdir = 'ion_scale/5x3x192x48x24/0T';
% resdir = 'ion_scale/5x3x192x48x24/no_gradN/0T';
% resdir = 'ion_scale/5x3x192x48x24/lower_grad/PT';
% resdir = 'ion_scale/9x5x256x64x32/0T';
% resdir = 'ion_scale/restart/5x3x256x64x32/0T';
% resdir = 'ion_scale/restart/9x5x192x48x24/0T';
% resdir = 'adiabatic_electrons/5x2x256x64x32/0T';
% resdir = 'adiabatic_electrons/5x2x192x48x24/0T';
% resdir = 'hot_electrons/256x64x32/0T';
% resdir = 'hot_electrons/256x64x32/0T';
% resdir = 'hot_electrons/512x64x32/0T';

% No grad N
% PARTITION = '/misc/gyacomo23_outputs/triangularity_paper/no_gradN';
% resdir = '/ion_scale/3x2x256x64x24/0T';
% resdir = '/ion_scale/5x2x256x64x24/0T';
% resdir = '/ion_scale/9x5x128x32x24/0T';
% resdir = '/adiabatic_electrons/3x2x256x64x24/0T';
% resdir = '/adiabatic_electrons/5x2x256x64x24/0T';
% resdir = '/hot_electrons/L_300/256x64x32/0T';

% PARTITION = '/home/ahoffman/gyacomo/results/lin_DIIID_LM_rho95_scan/';
% resdir   = '6x2x32_5x3_Lx_120_Ly_37.6991_q0_4.8_e_0.3_s_2.5_mil__kN_1.7_kT_5.2_nu_2.0e-02_DGGK/';
% resdir   = '6x2x32_17x9_Lx_120_Ly_37.6991_q0_4.8_e_0.3_s_2.5_mil__kN_1.7_kT_5.2_nu_2.0e-02_DGGK/';
% resdir   = '6x2x32_5x3_Lx_120_Ly_12.5664_q0_4.8_e_0.3_s_2.5_mil__kN_1.7_kT_5.2_nu_2.0e-02_DGGK/';
% resdir   = '6x2x32_17x9_Lx_120_Ly_12.5664_q0_4.8_e_0.3_s_2.5_mil__kN_1.7_kT_5.2_nu_2.0e-02_DGGK/';
% resdir   = '6x2x32_5x3_Lx_120_Ly_8.1955_q0_4.8_e_0.3_s_2.5_mil__kN_1.7_kT_5.2_nu_2.0e-02_DGGK/';
% resdir   = '6x2x32_17x9_Lx_120_Ly_8.1955_q0_4.8_e_0.3_s_2.5_mil__kN_1.7_kT_5.2_nu_2.0e-02_DGGK/';

DATADIR = [PARTITION,resdir,'/'];
% DATADIR = '/home/ahoffman/gyacomo/simulations/ralf/2D_Zpinch_ITG/3x2x64x48x1_no_curvB/';
% DATADIR = '/home/ahoffman/gyacomo/simulations/ralf/3D_Zpinch_ITG/3x2x64x48x16_nocurvB/';
% DATADIR = '/home/ahoffman/gyacomo/simulations/ralf/3D_Zpinch_ITG/3x2x64x48x16_nocurvB_-14/';
% DATADIR = '/home/ahoffman/gyacomo/simulations/ricci_UHD/';
read_flux_out_XX(DATADIR,1,1);
%%
J0 = 00; J1 = 10;

% Load basic info (grids and time traces)
data    = {};
data    = compile_results_low_mem(data,DATADIR,J0,J1);
[data.Na00, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'Na00');
data.Ni00 = reshape(data.Na00(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
try
data.Ne00 = reshape(data.Na00(2,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
catch
end
[data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');
if data.inputs.BETA > 0
    [data.PSI, data.Ts3D]  = compile_results_3D(DATADIR,J0,J1,'psi');
end
if 1
    %%
[data.TEMP, data.Ts3D] = compile_results_3Da(data.folder,J0,J1,'temp');
% [data.UPAR, data.Ts3D] = compile_results_3Da(data.folder,J0,J1,'upar');
% [data.UPER, data.Ts3D] = compile_results_3Da(data.folder,J0,J1,'uper');
[data.DENS, data.Ts3D] = compile_results_3Da(data.folder,J0,J1,'dens');
data.TEMP_I = reshape(data.TEMP(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
% data.UPAR_I = reshape(data.UPAR(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
% data.UPER_I = reshape(data.UPER(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
data.DENS_I = reshape(data.DENS(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
data.Ni00 = reshape(data.Na00(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
if data.inputs.Na > 1
    data.TEMP_E = reshape(data.TEMP(2,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
    data.DENS_E = reshape(data.DENS(2,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
    data.Ne00 = reshape(data.Na00(2,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
end
end
if 0
%% Plot transport and phi radial profile
% [data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');
% [data.PSI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'psi');
options.TAVG_0   = 250;
options.TAVG_1   = options.TAVG_0+50;
options.NCUT     = 5;              % Number of cuts for averaging and error estimation
options.NMVA     = 1;              % Moving average for time traces
% options.ST_FIELD = '\Gamma_x';   % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
% options.ST_FIELD = 'n_i';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
% options.ST_FIELD = 'u_i';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
% options.ST_FIELD = 'T_i';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
% options.ST_FIELD = 'n_i T_i';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.ST_FIELD = '\phi';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
% options.ST_FIELD = 'Q_{xi}';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
% options.ST_FIELD = 'G_x';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
% options.ST_FIELD = 'w_{Ez}';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
% options.ST_FIELD = 'v_{Ey}';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
% options.ST_FIELD = 'N_i^{00}';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.INTERP   = 0;
options.RESOLUTION = 256;
plot_radial_transport_and_spacetime(data,options);
end

if 1
%% 2D field snapshots
% Options
options.INTERP    = 0;
options.POLARPLOT = 0;
options.AXISEQUAL = 0;
options.NORMALIZE = 0;
options.LOGSCALE  = 0;
options.CLIMAUTO  = 1;
options.TAVG      = 1;
options.NAME      = ['N_i^{00}'];
% options.NAME      = 'n_i';
% options.NAME      = 'upar_i';
% options.NAME      = 'T_i';
% options.NAME      = 'Q_{xi}';
% options.NAME      = 'v_{Ey}';
% options.NAME      = 'w_{Ez}';
% options.NAME      = '\omega_z';
% options.NAME      = '\phi';
% options.NAME      = 'n_i-n_e';
loc =11;
[~,i_] = min(abs(loc - data.grids.y));
options.COMP =i_;
% options.PLAN      = '3D';  
options.PLAN      = 'xy'; options.COMP =floor(data.grids.Nz/2)+1; 
% options.PLAN      = 'xz'; options.COMP ='avg';
% options.COMP ='avg'; 
options.XYZ  =[-11 20 -2]; 
options.TIME = [0.0 0.1 0.2]; options.TAVG = 0;
% options.TIME = [50:500]; options.TAVG = 1;
options.RESOLUTION = 256;
fig = photomaton(data,options);
colormap(gray)
clim('auto')
% set(gca,'ColorScale','log')
% save_figure(data,fig)
end
if 0
%% Performance profiler
profiler(data)
end

if 0
%% Mode evolution
% [data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');
% [data.Na00, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'Na00');
% data.Ni00 = reshape(data.Na00(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
% data.Ne00 = reshape(data.Na00(2,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));

options.NORMALIZED = 0;
options.TIME   = data.Ts3D;
options.KX_TW  = [0.5 1]*data.Ts3D(end); %kx Growth rate time window
options.KY_TW  = [0.5 1]*data.Ts3D(end);  %ky Growth rate time window
options.NMA    = 1;
options.NMODES = 64;
options.iz     = 'avg'; % avg or index
options.ik     = 1; % sum, max or index
options.fftz.flag = 0;
options.FIELD  = 'Ni00';
% options.FIELD  = 'phi';
% options.FIELD  = 'T_i';
options.GOK2   = 0;
options.SHOWFIG = 1;
[fig, wkykx, ekykx] = mode_growth_meter(data,options);
% %%
% kx = (1:data.grids.Nx/2)'*2*pi/data.fort_00.GRID.Lx;
% ky = (1:data.grids.Ny/2)'*2*pi/data.fort_00.GRID.Ly;
% gkxky = real(wkykx(2:end,1:data.grids.Nx/2))';
% gkxky(isnan(gkxky)) =0;
% gkxky(isinf(gkxky)) =0;
% % gkxky(gkxky<0)      =0;
% % gkxky = imgaussfilt(gkxky,1);
% %
% wkxky = imag(wkykx(2:end,1:data.grids.Nx/2))';
% wkxky(isnan(wkxky)) =0;
% wkxky(isinf(wkxky)) =0;
% % wkxky(wkxky<0)      =0;
% % wkxky = imgaussfilt(wkxky,1.5);
% %
% figure; 
% subplot(121)
%     contourf(kx,ky,gkxky',10)
%     % clim(0.5*[0 1]); 
%     % colormap(bluewhitered); colorbar;
%     xlim([0.025 1]);
%     xlabel('$k_x\rho_s$'); ylabel('$k_y\rho_s$')
% subplot(122)
%     contourf(kx,ky,wkxky',10)
%     % clim(1*[0 1]); 
%     % colormap(bluewhitered); colorbar 
%     xlim([0.025 1]);
%     xlabel('$k_x\rho_s$'); ylabel('$k_y\rho_s$')
% % save_figure(data,fig,'.png')
end

if 0
%% Hermite-Laguerre spectrum
[data.Napjz, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'Napjz');
% data.Napjz(1,3,1,:,:) = data.Napjz(1,3,1,:,:)*data.inputs.tau(1);
% data.Napjz(1,1,2,:,:) = data.Napjz(1,1,2,:,:)*data.inputs.tau(1);
% [data.Napjz, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'Nipjz');
options.ST         = 1;
options.NORMALIZED = 0;
options.LOGSCALE   = 1;
options.FILTER     = 0; %filter the 50% time-average of the spectrum from
options.TAVG_2D    = 0; %Show a 2D plot of the modes, 50% time averaged
options.TAVG_2D_CTR= 0; %make it contour plot
options.TWINDOW    = [0 20];
fig = show_moments_spectrum(data,options);
end

if (0 && NZ>4)
%% Ballooning plot
% [data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');
if data.inputs.BETA > 0
[data.PSI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'psi');
end
options.time_2_plot = [25 100];
options.kymodes     = 1.5;
options.normalized  = 1;
options.PLOT_KP     = 0;
% options.field       = 'phi';
options.SHOWFIG  = 1;
[fig, chi, phib, psib, ~] = plot_ballooning(data,options);
end

if 0
%% 1D spectral plot
options.TIME  = [100 300]; % averaging time window
options.NAME      = ['N_i^{00}'];
% options.NAME      = 'n_i';
% options.NAME      = 'T_i';
% options.NAME      = 'Q_{xi}';
% options.NAME      = 's_{Ey}';
% options.NAME      = '\phi';
% options.NAME      = '\psi';
options.NORMALIZE = 0;
[fig] = plot_spectrum(data,options);
end

if 0
%% 1D radial plot
options.TIME  = [100 300]; % averaging time window
options.NAME      = ['N_i^{00}'];
% options.NAME      = 'n_i';
% options.NAME      = 'T_i';
% options.NAME      = 'Q_{xi}';
% options.NAME      = 's_{Ey}';
% options.NAME      = '\phi';
% options.NAME      = '\psi';
options.NORMALIZE = 0;
[fig] = plot_spectrum(data,options);
end


if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
% [data.PHI, data.Ts3D]  = compile_results_3D(DATADIR,J0,J1,'phi');
% [data.PSI, data.Ts3D]  = compile_results_3D(DATADIR,J0,J1,'psi');
% [data.Na00, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'Na00');
% data.Ni00 = reshape(data.Na00(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
% data.Ne00 = reshape(data.Na00(2,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
% [data.DENS, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'dens');
% data.DENS_I = reshape(data.DENS(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
% [data.TEMP, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'temp');
% data.TEMP_I = reshape(data.TEMP(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
options.INTERP    = 1;
options.POLARPLOT = 0;
options.BWR       = 0; % bluewhitered plot or gray
options.CLIMAUTO  = 1; % adjust the colormap auto
% options.NAME      = '\phi';
% options.NAME      = 'w_{Ez}';
% options.NAME      = '\psi';
% options.NAME      = 'T_i';
% options.NAME      = '\phi^{NZ}';
options.NAME     = ['N_i^{00}'];
% options.NAME     = ['N_e^{00}'];
options.PLAN      = 'xy'; options.COMP =floor(data.grids.Nz/2)+1; 
% options.PLAN      = 'xz'; options.COMP ='avg';
% options.PLAN      = '3D';  
options.XYZ  =[-21 20 0]; 
options.TIME      =  data.Ts3D(1:1:end);
% options.TIME      = [0:1500];
data.EPS          = 0.1;
data.a = data.EPS * 2000;
options.RESOLUTION = 256;
options.FPS       = 12;
options.RMAXIS    = 1;
create_film(data,options,'.gif');
end

if 0
%% Metric infos
options.SHOW_FLUXSURF = 1;
options.SHOW_METRICS  = 1;
[fig, geo_arrays] = plot_metric(data,options);
end

if 0
%% Study singular values
[data.SV_ky_pj, data.Ts2D] = compile_results_2D(DATADIR,J0,J1,'sv_ky_pj');
nSV = data.grids.Np * data.grids.Nj;
colors_ = jet(nSV);
figure
for i = 1:nSV
    sv = squeeze(data.SV_ky_pj(i,:));
    semilogy(data.Ts2D,sv,...
        'color',colors_(i,:),'DisplayName',['SV ',num2str(i)]);hold on
end
legend('show');
end

if 0
%% Pseudo fluid analysis
Time_window = [100 200];
pseudo_fluid_analysis
end