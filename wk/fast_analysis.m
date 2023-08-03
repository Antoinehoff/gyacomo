gyacomodir = pwd; gyacomodir = gyacomodir(1:end-2); % get code directory
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add
default_plots_options
% Partition of the computer where the data have to be searched
PARTITION='';

%% Paper 3
resdir = '/Users/ahoffmann/gyacomxlio/results/paper_3/DTT_rho85/3x2x192x48x32';
% resdir = '/Users/ahoffmann/gyacomo/results/paper_3/DTT_rho85/3x2x192x48x32_NT';
% resdir = '/Users/ahoffmann/gyacomo/results/paper_3/DTT_rho98/3x2x192x48x32';
resdir = '/Users/ahoffmann/gyacomo/results/paper_3/DTT_rho98/3x2x192x48x32_0.25grad';
 %%
J0 = 00; J1 = 02;

% Load basic info (grids and time traces)
DATADIR = [PARTITION,resdir,'/'];
data    = {};
data    = compile_results_low_mem(data,DATADIR,J0,J1);

if 1
%% Plot transport and phi radial profile
[data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');
[data.PSI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'psi');
options.TAVG_0   = 100;
options.TAVG_1   = 1000;
options.NCUT     = 5;              % Number of cuts for averaging and error estimation
options.NMVA     = 1;              % Moving average for time traces
% options.ST_FIELD = '\Gamma_x';   % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.ST_FIELD = '\phi';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.INTERP   = 0;
options.RESOLUTION = 256;
plot_radial_transport_and_spacetime(data,options);
end

if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
% [data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');
[data.Na00, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'Na00');
data.Ni00 = reshape(data.Na00(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
options.INTERP    = 1;
options.POLARPLOT = 0;
% options.NAME      = '\phi';
% options.NAME      = '\phi^{NZ}';
% options.NAME      = '\omega_z';
options.NAME     = 'N_i^{00}';
% options.NAME      = 's_{Ey}';
% options.NAME      = 'n_i^{NZ}';
% options.NAME      = 'Q_x';
% options.NAME      = 'n_i';
% options.NAME      = 'n_i-n_e';
options.PLAN      = 'kxky';
% options.NAME      = 'f_i';
% options.PLAN      = 'sx';
options.COMP      = 1;
% options.TIME      = data.Ts5D(end-30:end);
options.TIME      =  data.Ts3D(1:1:end);
% options.TIME      = [0:1500];
data.EPS          = 0.1;
data.a = data.EPS * 2000;
options.RESOLUTION = 256;
create_film(data,options,'.gif')
end

if 0
%% field snapshots
% Options
[data.Na00, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'Na00');
[data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');
data.Ni00 = reshape(data.Na00(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));

options.INTERP    = 0;
options.POLARPLOT = 0;
options.AXISEQUAL = 0;
options.NORMALIZE = 0;
% options.NAME      = 'N_i^{00}';
options.NAME      = '\phi';
options.PLAN      = 'kxky';
options.COMP      = 'avg';
options.TIME      = [0 10 20 50];
% options.TIME      = data.Ts3D(1:2:end);
options.RESOLUTION = 256;
fig = photomaton(data,options);
% save_figure(data,fig)
end
if 0
%% Performance profiler
profiler(data)
end

if 1
%% Hermite-Laguerre spectrum
[data.Napjz, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'Napjz');
% [data.Napjz, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'Nipjz');
options.ST         = 0;
options.NORMALIZED = 0;
options.LOGSCALE   = 1;
options.FILTER     = 0; %filter the 50% time-average of the spectrum from
options.TAVG_2D    = 0; %Show a 2D plot of the modes, 50% time averaged
options.TAVG_2D_CTR= 0; %make it contour plot
fig = show_moments_spectrum(data,options);
end

if 0
%% Mode evolution
[data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');

options.NORMALIZED = 0;
options.TIME   = [000:9000];
options.KX_TW  = [30 40]; %kx Growth rate time window
options.KY_TW  = [10 20];  %ky Growth rate time window
options.NMA    = 1;
options.NMODES = 800;
options.iz     = 'avg'; % avg or index
options.ik     = 1; % sum, max or index
options.fftz.flag = 0;
fig = mode_growth_meter(data,options);
% save_figure(data,fig,'.png')
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