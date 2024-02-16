gyacomodir = '../../../';
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add
default_plots_options

J0 = 00; J1 = 00;

% Load basic info (grids and time traces)
DATADIR = [pwd,'/'];
data    = {};
data    = compile_results_low_mem(data,DATADIR,J0,J1);
[data.Na00, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'Na00');
data.Ni00 = reshape(data.Na00(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));

% field snapshots
options.INTERP    = 0;
options.POLARPLOT = 0;
options.AXISEQUAL = 1;
options.NORMALIZE = 0;
options.LOGSCALE  = 0;
options.CLIMAUTO  = 1;
options.NAME      = ['N_i^{00}'];
options.PLAN      = 'xy'; options.COMP =floor(data.grids.Nz/2)+1;
options.TIME      = [0 2.0 4.0];
options.RESOLUTION = 256;
options.BWR       = 0; % bluewhitered plot or gray
fig = photomaton(data,options);
colormap(gray)
clim('auto')
data.FIGDIR = DATADIR;
% save_figure(data,fig,'.png');


if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.INTERP    = 1;
options.POLARPLOT = 0;
options.BWR       = 0; % bluewhitered plot or gray
options.CLIMAUTO  = 1; % adjust the colormap auto
options.NAME     = ['N_i^{00}'];
% options.NAME     = ['N_i^{00}'];
options.PLAN      = ['xy'];
options.COMP      = 1;
options.TIME      =  data.Ts3D(1:1:end);
data.EPS          = 0.1;
data.a = data.EPS * 2000;
options.RESOLUTION = 64;
options.FPS       = 12;
options.RMAXIS    = 1;
create_film(data,options,'.gif')
end
