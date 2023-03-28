% Directory of the code "mypathtogyacomo/gyacomo/"
% Partition of the computer where the data have to be searched
PARTITION  = '/misc/gyacomo23_outputs/';
% PARTITION  = gyacomodir;

%% CBC 
% resdir = 'paper_2_GYAC23/CBC/7x4x192x96x32_nu_0.05_muxy_0.5_muz_0.2';
% resdir = 'paper_2_GYAC23/CBC/7x4x192x96x32_nu_0.05_muxy_1.0_muz_1.0';
% resdir = 'paper_2_GYAC23/CBC/7x4x192x96x32_nu_0.05_muxy_1.0_muz_2.0';
% resdir = 'paper_2_GYAC23/CBC/Full_NL_7x4x192x96x32_nu_0.05_muxy_1.0_muz_2.0';

%% tests
resdir = 'paper_2_GYAC23/precision_study/5x3x128x64x24';
% resdir = 'paper_2_GYAC23/precision_study/5x3x128x64x24_Lx_180';
%%
J0 = 00; J1 = 10;

% Load basic info (grids and time traces)
DATADIR = [PARTITION,resdir,'/'];
data    = {};
data    = compile_results_low_mem(data,DATADIR,J0,J1);

%% Plot transport and phi radial profile
[data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');

options.TAVG_0   = 100;
options.TAVG_1   = 1000;
options.NCUT     = 5;              % Number of cuts for averaging and error estimation
options.NMVA     = 1;              % Moving average for time traces
% options.ST_FIELD = '\Gamma_x';   % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.ST_FIELD = '\phi';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.INTERP   = 0;
options.RESOLUTION = 256;
fig = plot_radial_transport_and_spacetime(data,options);

if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
options.INTERP    = 0;
options.POLARPLOT = 0;
options.NAME      = '\phi';
% options.NAME      = '\omega_z';
% options.NAME     = 'N_i^{00}';
% options.NAME      = 's_{Ey}';
% options.NAME      = 'n_i^{NZ}';
% options.NAME      = 'Q_x';
% options.NAME      = 'n_i';
% options.NAME      = 'n_i-n_e';
options.PLAN      = 'xy';
% options.NAME      = 'f_i';
% options.PLAN      = 'sx';
options.COMP      = 'avg';
% options.TIME      = data.Ts5D(end-30:end);
% options.TIME      =  data.Ts3D;
options.TIME      = [0:1500];
data.EPS          = 0.1;
data.a = data.EPS * 2000;
options.RESOLUTION = 256;
create_film(data,options,'.gif')
end
