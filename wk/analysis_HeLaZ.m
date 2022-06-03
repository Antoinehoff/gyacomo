addpath(genpath([helazdir,'matlab'])) % ... add
addpath(genpath([helazdir,'matlab/plot'])) % ... add
addpath(genpath([helazdir,'matlab/compute'])) % ... add
addpath(genpath([helazdir,'matlab/load'])) % ... add

%% Load the results
LOCALDIR  = [helazdir,'results/',outfile,'/'];
MISCDIR   = ['/misc/HeLaZ_outputs/results/',outfile,'/'];
system(['mkdir -p ',MISCDIR]);
CMD = ['rsync ', LOCALDIR,'outputs* ',MISCDIR]; disp(CMD);
system(CMD);
% Load outputs from jobnummin up to jobnummax
JOBNUMMIN = 10; JOBNUMMAX = 20;
data = compile_results(MISCDIR,JOBNUMMIN,JOBNUMMAX); %Compile the results from first output found to JOBNUMMAX if existing
data.localdir = LOCALDIR;

%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default_plots_options
disp('Plots')
FMT = '.fig';

if 0
%% Space time diagramm (fig 11 Ivanov 2020)
options.TAVG_0   = 0.98*data.Ts3D(end);
options.TAVG_1   = data.Ts3D(end); % Averaging times duration
options.NMVA     = 1;              % Moving average for time traces
% options.ST_FIELD = '\Gamma_x';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.ST_FIELD = '\phi';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.INTERP   = 1;
fig = plot_radial_transport_and_spacetime(data,options);
save_figure(data,fig)
end

if 0
%% statistical transport averaging
options.T = [16000 17000];
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
options.PLAN      = 'xz';
% options.NAME      = 'f_e';
% options.PLAN      = 'sx';
options.COMP      = 'avg';
% options.TIME      = dat.Ts5D;
options.TIME      = 1250:1:1500;
data.EPS          = 0.1;
data.a = data.EPS * 2000;
create_film(data,options,'.gif')
end

if 0
%% 2D snapshots
% Options
options.INTERP    = 0;
options.POLARPLOT = 0;
options.AXISEQUAL = 1;
options.NAME      = '\phi';
% options.NAME      = 'n_i';
% options.NAME      = 'N_i^{00}';
% options.NAME      = 'T_i';
% options.NAME      = '\Gamma_x';
% options.NAME      = 'k^2n_e';
options.PLAN      = 'kxky';
% options.NAME      = 'f_i';
% options.PLAN      = 'sx';
options.COMP      = 'avg';
options.TIME      = [1200 1300 1400 1500];
data.a = data.EPS * 2e3;
fig = photomaton(data,options);
save_figure(data,fig)
end

if 1
%% Ballooning plot
options.time_2_plot = [800 900];
options.kymodes     = [0.5];
options.normalized  = 1;
options.sheared     = 0;
options.field       = 'phi';
fig = plot_ballooning(data,options);
end

if 0
%% Kinetic distribution function sqrt(<f_a^2>xy) (GENE vsp)
% options.SPAR      = linspace(-3,3,64)+(6/127/2);
% options.XPERP     = linspace( 0,6,64);
options.SPAR      = gene_data.vp';
options.XPERP     = gene_data.mu';
options.iz        = 9;
options.T         = 30;
options.PLT_FCT   = 'contour';
options.ONED      = 0;
options.non_adiab = 1;
options.SPECIE    = 'i';
options.RMS       = 1; % Root mean square i.e. sqrt(sum_k|f_k|^2) as in Gene
fig = plot_fa(data,options);
save_figure(data,fig)
end

if 0
%% Hermite-Laguerre spectrum
% options.TIME = 'avg';
options.P2J        = 1;
options.ST         = 0;
options.PLOT_TYPE  = 'space-time';
options.NORMALIZED = 1;
options.JOBNUM     = 0;
options.TIME       = [1300 1500];
options.specie     = 'i';
options.compz      = 'avg';
fig = show_moments_spectrum(data,options);
% fig = show_napjz(data,options);
save_figure(data,fig)
end

if 0
%% Time averaged spectrum
options.TIME   = 1000:1200;
options.NORM   =1;
options.NAME   = '\phi';
% options.NAME      = 'n_i';
% options.NAME      ='\Gamma_x';
options.PLAN   = 'kxky';
options.COMPZ  = 'avg';
options.OK     = 0;
options.COMPXY = 'avg';
options.COMPT  = 'avg';
options.PLOT   = 'semilogy';
fig = spectrum_1D(data,options);
% save_figure(data,fig)
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
% save_figure(data,fig)
end

if 0
%% Mode evolution
options.NORMALIZED = 1;
options.K2PLOT = 1;
options.TIME   = 5:30;
options.NMA    = 1;
options.NMODES = 15;
options.iz     = 9;
fig = mode_growth_meter(data,options);
save_figure(gbms_dat,fig)
end

if 0
%% ZF caracteristics (space time diagrams)
TAVG_0 = 200; TAVG_1 = 3000; % Averaging times duration
% chose your field to plot in spacetime diag (uzf,szf,Gx)
fig = ZF_spacetime(data,TAVG_0,TAVG_1);
save_figure(data,fig)
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
save_figure(data,fig)
end

if 0
%% 3D plot on the geometry
options.INTERP    = 1;
options.NAME      = 'n_i';
options.PLANES    = [1:2:12];
options.TIME      = [60];
options.PLT_MTOPO = 1;
data.rho_o_R      = 2e-3; % Sound larmor radius over Machine size ratio
fig = show_geometry(data,options);
save_figure(data,fig)
end
