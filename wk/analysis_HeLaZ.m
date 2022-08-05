addpath(genpath([helazdir,'matlab'])) % ... add
addpath(genpath([helazdir,'matlab/plot'])) % ... add
addpath(genpath([helazdir,'matlab/compute'])) % ... add
addpath(genpath([helazdir,'matlab/load'])) % ... add

%% Load the results
LOCALDIR  = [helazdir,'results/',outfile,'/'];
MISCDIR   = ['/misc/HeLaZ_outputs/results/',outfile,'/'];
system(['mkdir -p ',MISCDIR]);
system(['mkdir -p ',LOCALDIR]);
CMD = ['rsync ', LOCALDIR,'outputs* ',MISCDIR]; disp(CMD);
system(CMD);
% Load outputs from jobnummin up to jobnummax
data = compile_results(MISCDIR,JOBNUMMIN,JOBNUMMAX); %Compile the results from first output found to JOBNUMMAX if existing
data.localdir = LOCALDIR;
data.FIGDIR   = LOCALDIR;

%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default_plots_options
disp('Plots')
FMT = '.fig';

if 1
%% Space time diagramm (fig 11 Ivanov 2020)
% data.scale = 1;%/(data.Nx*data.Ny)^2;
options.TAVG_0   = 25;%0.4*data.Ts3D(end); 
options.TAVG_1   = 40;%0.9*data.Ts3D(end); % Averaging times duration
options.NCUT     = 4;              % Number of cuts for averaging and error estimation
options.NMVA     = 1;              % Moving average for time traces
% options.ST_FIELD = '\Gamma_x';   % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.ST_FIELD = '\phi';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.INTERP   = 1;
fig = plot_radial_transport_and_spacetime(data,options);
save_figure(data,fig,'.png')
end

if 0
%% statistical transport averaging
options.T = [200 400];
fig = statistical_transport_averaging(data,options);
end
    
if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
options.INTERP    = 1;
options.POLARPLOT = 0;
options.NAME      = '\phi';
% options.NAME      = 'N_i^{00}';
% options.NAME      = 'v_y';
% options.NAME      = 'n_i^{NZ}';
% options.NAME      = '\Gamma_x';
% options.NAME      = 'n_i';
options.PLAN      = 'kxky';
% options.NAME      = 'f_i';
% options.PLAN      = 'sx';
options.COMP      = 'avg';
% options.TIME      = data.Ts5D(end-30:end);
options.TIME      =  data.Ts3D;
% options.TIME      = [850:0.1:1000];
data.EPS          = 0.1;
data.a = data.EPS * 2000;
create_film(data,options,'.gif')
end

if 0
%% 2D snapshots
% Options
options.INTERP    = 1;
options.POLARPLOT = 0;
options.AXISEQUAL = 0;
options.NAME      = '\phi';
% options.NAME      = 'n_e';
% options.NAME      = 'N_i^{00}';
% options.NAME      = 'T_i';
% options.NAME      = '\Gamma_x';
% options.NAME      = 'k^2n_e';
options.PLAN      = 'kxky';
% options.NAME      'f_i';
% options.PLAN      = 'sx';
options.COMP      = 'avg';
options.TIME      = [400 440];
data.a = data.EPS * 2e3;
fig = photomaton(data,options);
% save_figure(data,fig)
end

if 0
%% 3D plot on the geometry
options.INTERP    = 1;
options.NAME      = '\phi';
options.PLANES    = [1:1:16];
options.TIME      = [15];
options.PLT_MTOPO = 0;
data.rho_o_R      = 2e-3; % Sound larmor radius over Machine size ratio
fig = show_geometry(data,options);
save_figure(data,fig,'.png')
end

if 0
%% Kinetic distribution function sqrt(<f_a^2>xy) (GENE vsp)
% options.SPAR      = linspace(-3,3,32)+(6/127/2);
% options.XPERP     = linspace( 0,6,32);
options.SPAR      = gene_data.vp';
options.XPERP     = gene_data.mu';
options.iz        = 'avg';
options.T         = [300 600];
options.PLT_FCT   = 'contour';
options.ONED      = 0;
options.non_adiab = 0;
options.SPECIE    = 'i';
options.RMS       = 1; % Root mean square i.e. sqrt(sum_k|f_k|^2) as in Gene
fig = plot_fa(data,options);
% save_figure(data,fig,'.png')
end

if 0
%% Hermite-Laguerre spectrum
% options.TIME = 'avg';
options.P2J        = 0;
options.ST         = 1;
options.PLOT_TYPE  = 'space-time';
options.NORMALIZED = 0;
options.JOBNUM     = 0;
options.TIME       = [1000];
options.specie     = 'i';
options.compz      = 'avg';
fig = show_moments_spectrum(data,options);
% fig = show_napjz(data,options);
save_figure(data,fig,'.png');
end

if 0
%% Time averaged spectrum
options.TIME   = [300 600];
options.NORM   =1;
options.NAME   = '\phi';
% options.NAME      = 'N_i^{00}';
% options.NAME   ='\Gamma_x';
options.PLAN   = 'kxky';
options.COMPZ  = 'avg';
options.OK     = 0;
options.COMPXY = 'avg'; % avg/sum/max/zero/ 2D plot otherwise
options.COMPT  = 'avg';
options.PLOT   = 'semilogy';
fig = spectrum_1D(data,options);
% save_figure(data,fig,'.png')
end

if 0
%% 1D real plot
options.TIME   = [50 100 200];
options.NORM   = 0;
options.NAME   = '\phi';
% options.NAME      = 'n_i';
% options.NAME      ='\Gamma_x';
% options.NAME      ='s_y';
options.COMPX  = 'avg';
options.COMPY  = 'avg';
options.COMPZ  = 1;
options.COMPT  = 1;
options.MOVMT  = 1;
fig = real_plot_1D(data,options);
% save_figure(data,fig,'.png')
end

if 0
%% Mode evolution
options.NORMALIZED = 0;
options.K2PLOT = 1;
options.TIME   = [0:160];
options.NMA    = 1;
options.NMODES = 1;
options.iz     = 'avg';
fig = mode_growth_meter(data,options);
save_figure(data,fig,'.png')
end

if 0
%% ZF caracteristics (space time diagrams)
TAVG_0 = 1200; TAVG_1 = 1500; % Averaging times duration
% chose your field to plot in spacetime diag (uzf,szf,Gx)
fig = ZF_spacetime(data,TAVG_0,TAVG_1);
save_figure(data,fig,'.png')
end

if 0
%% Metric infos
fig = plot_metric(data);
end

if 0
%% linear growth rate for 3D fluxtube
trange = [0 100];
nplots = 1;
lg = compute_fluxtube_growth_rate(data,trange,nplots);
end

if 0
%% linear growth rate for 3D Zpinch
trange = [5 15];
options.keq0 = 1; % chose to plot planes at k=0 or max
options.kxky = 1;
options.kzkx = 0;
options.kzky = 1;
[lg, fig] = compute_3D_zpinch_growth_rate(data,trange,options);
save_figure(data,fig,'.png')
end
