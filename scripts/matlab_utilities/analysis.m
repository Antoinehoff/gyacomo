%% This is a example of matlab analysis script
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Before anything, setup the right path to your gyacommo directory here
gyacomodir = '/Users/ahoffmann/gyacomo/';
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% add to our path the routines in our local matlab folder
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add
% setup font and line size
default_plots_options


% -------------------------------------------------------------------------
% Indicate where are the data you would like to analyze
DATADIR = '/Users/ahoffmann/gyacomo/simulations/problem_01/';
% the rest is optional
% -------------------------------------------------------------------------

% Jobs to load (from outputs_J0.h5 to outputs_J1.h5 if exists)
J0 = 00; J1 = 10;

% main structure that will contain all the loaded data
data    = {};

% Load basic info (grids and time traces)
data    = compile_results_low_mem(data,DATADIR,J0,J1);

% load EM fields ----------------------------------------------------------
% Electrostatic potential
[data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');
% Electromagnetic potential (if beta non zero)
if data.inputs.BETA > 0
    [data.PSI, data.Ts3D]  = compile_results_3D(DATADIR,J0,J1,'psi');
end

% load moments ------------------------------------------------------------
% temperature
[TEMP, data.Ts3D] = compile_results_3Da(data.folder,J0,J1,'temp');
% parallel velocity
[UPAR, data.Ts3D] = compile_results_3Da(data.folder,J0,J1,'upar');
% perpendicular velocity
[UPER, data.Ts3D] = compile_results_3Da(data.folder,J0,J1,'uper');
% density
[DENS, data.Ts3D] = compile_results_3Da(data.folder,J0,J1,'dens');
% gyrocenter density (first GM)
[Na00, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'Na00');
% asign them to species
% first index is ions
data.TEMP_I = reshape(TEMP(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
data.UPAR_I = reshape(UPAR(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
data.UPER_I = reshape(UPER(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
data.DENS_I = reshape(DENS(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
data.Ni00   = reshape(Na00(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
% if we have a second species we store in electrons
if data.inputs.Na > 1
    data.TEMP_E = reshape(TEMP(2,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
    data.DENS_E = reshape(DENS(2,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
    data.Ne00 = reshape(Na00(2,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));
end
clear TEMP UPAR UPER DENS Na00;

%% Plot transport and phi radial profile
options.TAVG_0   = 25;              % averaging start time for the flux
options.TAVG_1   = 50;              % end time
options.NCUT     = 5;               % Number of cuts for averaging and error estimation
options.NMVA     = 1;               % Moving average for time traces
options.ST_FIELD = '\phi';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x,upar,N_i^{00} etc.)
options.INTERP   = 0;               % interp the pcolor plot or not
plot_radial_transport_and_spacetime(data,options);

%% The rest is optional
if 0
%% 2D field snapshots
% Options
options.INTERP    = 1;
options.AXISEQUAL = 0;
options.NORMALIZE = 0;
options.LOGSCALE  = 0;
% chose the name (available : n_i,upar_i,T_i,Q_{xi},v_{Ey},w_{Ez},\phi,\psi
options.NAME      = 'n_i';
% options.PLAN      = 'xy'; options.COMP =floor(data.grids.Nz/2)+1; 
options.PLAN      = 'kxky'; options.COMP =floor(data.grids.Nz/2)+1; 
% options.PLAN      = 'xz'; options.COMP ='avg';
% options.PLAN      = '3D';
options.XYZ  =[-11 20 -2]; 
options.TIME = [5 10 25 75]; options.TAVG = 0;
% options.TIME = [50:500]; options.TAVG = 1;
options.RESOLUTION = 256;
fig = photomaton(data,options);
colormap(gray)
% colorbar
% set(gca,'ColorScale','log')
% save_figure(data,fig)
end

if 0
%% Mode evolution
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
options.GOK2   = 0;
options.SHOWFIG = 1;
[fig, wkykx, ekykx] = mode_growth_meter(data,options);
end

if 0
%% Hermite-Laguerre spectrum
[data.Napjz, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'Napjz');
options.ST         = 1;
options.NORMALIZED = 0;
options.LOGSCALE   = 1;
options.FILTER     = 0; %filter the 50% time-average of the spectrum from
options.TAVG_2D    = 0; %Show a 2D plot of the modes, 50% time averaged
options.TAVG_2D_CTR= 0; %make it contour plot
options.TWINDOW    = [20 20];
fig = show_moments_spectrum(data,options);
end

if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.INTERP    = 1;
options.POLARPLOT = 0;
options.BWR       = 0; % bluewhitered plot or gray
options.CLIMAUTO  = 1; % adjust the colormap auto
options.NAME      = '\phi';
% options.NAME      = 'w_{Ez}';
% options.NAME      = '\psi';
% options.NAME      = 'T_i';
% options.NAME      = '\phi^{NZ}';
% options.NAME     = ['N_i^{00}'];
% options.NAME     = ['N_e^{00}'];
options.PLAN      = 'xy'; options.COMP =floor(data.grids.Nz/2)+1; 
% options.PLAN      = 'xz'; options.COMP ='avg';
% options.PLAN      = '3D';  
options.XYZ  =[-21 20 0]; 
options.TIME      =  data.Ts3D(1:1:end);
data.EPS          = 0.1;
data.a = data.EPS * 2000;
options.RESOLUTION = 256;
options.FPS       = 12;
options.RMAXIS    = 0;
create_film(data,options,'.gif');
end
