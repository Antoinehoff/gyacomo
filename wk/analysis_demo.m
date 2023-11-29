gyacomodir = pwd; gyacomodir = gyacomodir(1:end-2); % get code directory
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add
default_plots_options
options.RESOLUTION = 256;
% Partition of the computer where the data have to be searched

J0 = 00; J1 = 10;

% DATADIR = '/Users/ahoffmann/gyacomo/testcases/DIII-D_triangularity_fast_nonlinear/';
% DATADIR = '/Users/ahoffmann/gyacomo/testcases/HEL_DIII-D_triangularity/';
DATADIR = '/Users/ahoffmann/gyacomo/testcases/cyclone_example/';
data    = {};
data    = compile_results_low_mem(data,DATADIR,J0,J1);

if 1
%% Plot transport and phi radial profile
[data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');
options.TAVG_0   = 100;
options.TAVG_1   = 1000;
options.NCUT     = 5;              % Number of cuts for averaging and error estimation
options.NMVA     = 1;              % Moving average for time traces
options.ST_FIELD = '\phi';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.INTERP   = 1;
plot_radial_transport_and_spacetime(data,options);
end

if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
[data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');
[data.Na00, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'Na00');
data.Ni00 = reshape(data.Na00(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
options.INTERP    = 1;
options.POLARPLOT = 0;
options.BWR       = 1; % bluewhitered plot or gray
options.CLIMAUTO  = 0; % adjust the colormap auto
options.NAME      = '\phi';
% options.NAME     = 'N_i^{00}';
options.PLAN      = 'xy';
options.COMP      = 9;
options.TIME      =  data.Ts3D(1:2:end);
data.EPS          = 0.1;
data.a = data.EPS * 2000;
create_film(data,options,'.gif')
end

if 0
%% field snapshots
% Options
[data.Na00, data.Ts3D] = compile_results_3Da(data.folder,J0,J1,'Na00');
[data.PHI, data.Ts3D] = compile_results_3D(data.folder,J0,J1,'phi');
data.Ni00 = reshape(data.Na00(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));

options.INTERP    = 0;
options.POLARPLOT = 0;
options.AXISEQUAL = 0;
options.NORMALIZE = 0;
options.LOGSCALE  = 0;
options.CLIMAUTO  = 1;
options.NAME      = 'N_i^{00}';
% options.NAME      = 's_{Ey}';
% options.NAME      = '\phi';
options.PLAN      = 'yz';
options.COMP      = 'avg';
options.TIME      = [10 30];
fig = photomaton(data,options);
colormap(gray)
clim('auto')
% save_figure(data,fig)
end
if 0
%% Performance profiler
profiler(data)
end

if 0
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
