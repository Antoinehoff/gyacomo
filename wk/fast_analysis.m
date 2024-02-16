gyacomodir = pwd; gyacomodir = gyacomodir(1:end-2); % get code directory
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add
default_plots_options
% Partition of the computer where the data have to be searched
% PARTITION='/Users/ahoffmann/gyacomo/results/paper_3/';
PARTITION='/misc/gyacomo23_outputs/paper_3/';
% PARTITION = '../results/paper_3/';
%% Paper 3
% resdir = 'DTT_rho85/3x2x192x48x32';
% resdir = 'DTT_rho85/3x2x192x48x32_NT';
% resdir = 'DTT_rho98/3x2x192x48x32';
% resdir = 'DTT_rho98/3x2x192x48x32_0.25grad';
% resdir = 'LM_DIIID_rho95/5x3x512x92x32';
% resdir = 'LM_DIIID_rho95/3x2x512x92x32';
% resdir = 'DIIID_LM_rho90/3x2x256x128x32';
% resdir = 'DTT_rho85_geom_scan/P8_J4_delta_nuDGGK_conv_test/delta_-0.3_nu_0.9';
% resdir = 'NT_DIIID_Austin2019_rho95/3x2x256x64x32';

% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_2500/largerbox_moremodes';
% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_2500/PT/3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_2500/NT/3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_2500/VNT/3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_2500/VPT/3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_2500/CIRC/3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_2500/ELONG/3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_2500/PT/lin_3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_2500/NT/lin_3x2x128x32x32';

% resdir = 'DIIID_rho95_cold_ions_tau1e-3/huge';

% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_1000/PT/3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_1000/NT/3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_1000/PT/lin_3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_1000/NT/lin_3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_750/PT/lin_3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_750/NT/lin_3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_250/PT/lin_3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e-3/RT_250/NT/lin_3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e-3/huge';

% resdir = 'DIIID_rho95_cold_ions_tau1e-2/PT/3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e-2/NT/3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e0/PT/3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e0/PT/5x3x192x48x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e0/NT/3x2x128x32x32';
% resdir = 'DIIID_rho95_cold_ions_tau1e0/NT/5x3x192x48x32';
% resdir = 'DTT_rho85_cold_ions_tau1e-3/NT/3x2x128x32x32';
% % resdir = '../testcases/cyclone_example';
% resdir = '../testcases/CBC_ExBshear/std';
% resdir = '../results/paper_3/HM_DTT_rho98/3x2x128x64x64';

% resdir = 'DIIID_cold_ions_rho95_geom_scan/3x2x192x48x32_RT_1000_eps_q0_scan/NT/eps_0.35_q0_4.0';
% resdir = 'DTT_rho85_geom_scan/P2_J1_PT_sfact_shear_scan/shear_2.7_q0_-2.9';
% PARTITION = ''; resdir ='../results/HEL_CBC/tau_1e-3_kT_3500/128x32x24';
% PARTITION = ''; resdir ='../results/HEL_CBC/tau_1e-3_kT_3500/256x64x24';
% PARTITION = ''; resdir ='../results/HEL_CBC/tau_1e-3_colless';
% PARTITION = ''; resdir ='../results/HEL_CBC/tau_1e-3_kT_2000/lin_128x32x24';
% PARTITION = ''; resdir ='../results/HEL_CBC/tau_1e-3_kT_2000/128x32x24';
% PARTITION = ''; resdir ='../results/HEL_CBC/CBC_21/128x32x24';
% PARTITION = ''; resdir ='../results/HEL_CBC/192x48x24';
% PARTITION = ''; resdir ='../testcases/Ivanov_2020';
% PARTITION = ''; resdir ='../testcases/Hasegawa_Wakatani';
% resdir = 'HEL_CBC/256x92x24_max_trunc';

% resdir ='HEL_CBC/192x48x24';
% resdir ='HEL_CBC/256x92x24'3;
% resdir ='HEL_CBC/256x256x32/k_N__0.0_k_T__1750';

PARTITION = '/misc/gyacomo23_outputs/paper_3/DIIID_rho_95_Tstudy/';
% resdir = 'multi_scale_3x2x512x128x24';
% resdir = 'multi_scale_3x2x512x128x24_larger_box';
% resdir = 'multi_scale_3x2x768x192x24/continue_with_gradN_and_tau';
% resdir = 'multi_scale/multi_scale_5x2x768x192x24/PT';
% resdir = 'NT';
% resdir = 'electron_scale_3x2x256x64x24';
% resdir = 'electron_scale_3x2x128x64x24';
% resdir = 'ion_scale_3x2x192x48x32_larger_box';
% resdir = 'ion_scale_3x2x256x64x32'; % PT and NT also
% resdir = 'ion_scale/ion_scale_5x2x256x64x32_tau_1_RN_0/NT';
% resdir = 'adiab_e/5x2x256x64x32_tau_1_RN_0/NT';

resdir = 'multi_scale/multi_scale_3x2x768x192x24/NT';
% resdir = 'multi_scale/multi_scale_3x2x512x128x24_larger_box/NT';
% resdir = 'multi_scale/multi_scale_3x2x768x192x24/continue_with_gradN_and_tau';
% resdir = 'ion_scale/ion_scale_5x2x256x64x32_tau_1_RN_0/NT';
% resdir = 'adiab_e/5x2x256x64x32_tau_1_RN_0/0T';
% resdir = 'hot_electrons/hot_electrons_256x64x24/0T';



% PARTITION = '/misc/gyacomo23_outputs/paper_3/';
% resdir = 'DIIID_HEL_rho95/PT';

DATADIR = [PARTITION,resdir,'/'];
%%
J0 = 00; J1 = 20;

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
if 1
    %%
[data.TEMP, data.Ts3D] = compile_results_3Da(data.folder,J0,J1,'temp');
% [data.UPAR, data.Ts3D] = compile_results_3Da(data.folder,J0,J1,'upar');
[data.DENS, data.Ts3D] = compile_results_3Da(data.folder,J0,J1,'dens');
data.TEMP_I = reshape(data.TEMP(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
% data.UPAR_I = reshape(data.UPAR(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
data.DENS_I = reshape(data.DENS(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
data.Ni00 = reshape(data.Na00(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
if data.inputs.Na > 1
    data.TEMP_E = reshape(data.TEMP(2,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
    data.DENS_E = reshape(data.DENS(2,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
    data.Ne00 = reshape(data.Na00(2,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
end
end
if 1
%% Plot transport and phi radial profile
% [data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');
% [data.PSI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'psi');
options.TAVG_0   = data.Ts3D(end)/2;
options.TAVG_1   = data.Ts3D(end);
options.NCUT     = 5;              % Number of cuts for averaging and error estimation
options.NMVA     = 1;              % Moving average for time traces
% options.ST_FIELD = '\Gamma_x';   % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
% options.ST_FIELD = 'n_e';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
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
% options.NAME      = ['N_e^{00}'];
% options.NAME      = 'n_e';
% % options.NAME      = 'u_i';
% options.NAME      = 'T_i';
% options.NAME      = 'Q_{xe}';
options.NAME      = 'v_{Ey}';
% options.NAME      = 'w_{Ez}';
% options.NAME      = '\omega_z';
% options.NAME      = '\phi';
% options.NAME      = 'n_i-n_e';
loc =11;
[~,i_] = min(abs(loc - data.grids.y));
options.COMP =i_;
% options.PLAN      = '3D';  
options.PLAN      = 'xy'; options.COMP =floor(data.grids.Nz/2)+1; 
% options.PLAN      = 'yz'; options.COMP ='avg';
% options.COMP ='avg'; 
options.XYZ  =[-11 20 0]; 
options.TIME      = [90:150];
options.RESOLUTION = 256;
fig = photomaton(data,options);
% colormap(gray)
clim('auto')
% save_figure(data,fig)
end
if 0
%% Performance profiler
profiler(data)
end

if 1
%% Mode evolution
% [data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');
% [data.Na00, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'Na00');
% data.Ni00 = reshape(data.Na00(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
% data.Ne00 = reshape(data.Na00(2,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));

options.NORMALIZED = 0;
options.TIME   = data.Ts3D;
options.KX_TW  = [0.1 2.5]; %kx Growth rate time window
options.KY_TW  = [0.1 2.5];  %ky Growth rate time window
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
%%
kx = (1:data.grids.Nx/2)'*2*pi/data.fort_00.GRID.Lx;
ky = (1:data.grids.Ny/2)'*2*pi/data.fort_00.GRID.Ly;
gkxky = real(wkykx(2:end,1:data.grids.Nx/2))';
gkxky(isnan(gkxky)) =0;
gkxky(isinf(gkxky)) =0;
% gkxky(gkxky<0)      =0;
% gkxky = imgaussfilt(gkxky,1);
%
wkxky = imag(wkykx(2:end,1:data.grids.Nx/2))';
wkxky(isnan(wkxky)) =0;
wkxky(isinf(wkxky)) =0;
% wkxky(wkxky<0)      =0;
% wkxky = imgaussfilt(wkxky,1.5);
%
figure; 
subplot(121)
    contourf(kx,ky,gkxky',10)
    % clim(0.5*[0 1]); 
    % colormap(bluewhitered); colorbar;
    xlim([0.025 1]);
    xlabel('$k_x\rho_s$'); ylabel('$k_y\rho_s$')
subplot(122)
    contourf(kx,ky,wkxky',10)
    % clim(1*[0 1]); 
    % colormap(bluewhitered); colorbar 
    xlim([0.025 1]);
    xlabel('$k_x\rho_s$'); ylabel('$k_y\rho_s$')
% save_figure(data,fig,'.png')
end

if 1
%% Hermite-Laguerre spectrum
[data.Napjz, data.Ts3D] = compile_results_3Da(DATADIR,0,10,'Napjz');
% [data.Napjz, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'Nipjz');
options.ST         = 1;
options.NORMALIZED = 0;
options.LOGSCALE   = 0;
options.FILTER     = 0; %filter the 50% time-average of the spectrum from
options.TAVG_2D    = 0; %Show a 2D plot of the modes, 50% time averaged
options.TAVG_2D_CTR= 0; %make it contour plot
options.TWINDOW    = [6 20];
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

if 1
%% 1D spectral plot
options.TIME  = [30 80]; % averaging time window
% options.NAME      = ['N_i^{00}'];
% options.NAME      = 'n_i';
% options.NAME      = 'T_i';
% options.NAME      = 'Q_{xi}';
% options.NAME      = 's_{Ey}';
options.NAME      = '\psi';
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
options.INTERP    = 0;
options.POLARPLOT = 0;
options.BWR       = 1; % bluewhitered plot or gray
options.CLIMAUTO  = 0; % adjust the colormap auto
% options.NAME      = '\phi';
% options.NAME      = 'w_{Ez}';
% options.NAME      = '\psi';
options.NAME      = 'T_i';
% options.NAME      = '\phi^{NZ}';
% options.NAME     = ['N_e^{00}'];
% options.NAME     = ['N_i^{00}'];
options.PLAN      = 'xy';
% options.PLAN      = '3D';  
% options.XYZ  =[-11 20 0]; 
% options.COMP      = 'avg';
options.COMP      = floor(data.grids.Nz/2+1);
options.TIME      =  data.Ts3D(1:1:end);
% options.TIME      = [0:1500];
data.EPS          = 0.1;
data.a = data.EPS * 2000;
options.RESOLUTION = 256;
options.FPS       = 12;
options.RMAXIS    = 0;
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