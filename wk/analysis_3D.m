addpath(genpath('../matlab')) % ... add
addpath(genpath('../matlab/plot')) % ... add
addpath(genpath('../matlab/compute')) % ... add
addpath(genpath('../matlab/load')) % ... add

%% Load the results
LOCALDIR  = ['../results/',outfile,'/'];
MISCDIR   = ['/misc/HeLaZ_outputs/results/',outfile,'/'];
system(['mkdir -p ',MISCDIR]);
CMD = ['rsync ', LOCALDIR,'outputs* ',MISCDIR]; disp(CMD);
system(CMD);
% Load outputs from jobnummin up to jobnummax
JOBNUMMIN = 00; JOBNUMMAX = 20; 
data = compile_results(MISCDIR,JOBNUMMIN,JOBNUMMAX); %Compile the results from first output found to JOBNUMMAX if existing


%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default_plots_options
disp('Plots')
FMT = '.fig';

if 1
%% Space time diagramm (fig 11 Ivanov 2020)
TAVG_0 = 1000; TAVG_1 = 11000; % Averaging times duration
% chose your field to plot in spacetime diag (uzf,szf,Gx)
fig = plot_radial_transport_and_spacetime(data,TAVG_0,TAVG_1,'v_y',1);
save_figure(data,fig)
end

if 0
%% statistical transport averaging
options.T = [1500 2500];
fig = statistical_transport_averaging(data,options);
end

if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
options.INTERP    = 1;
options.POLARPLOT = 0;
options.NAME      = '\phi';
% options.NAME      = 'v_x';
% options.NAME      = 'n_i^{NZ}';
% options.NAME      = '\Gamma_x';
% options.NAME      = 'n_i';
% options.PLAN      = 'xy';
% options.NAME      = 'f_e';
% options.PLAN      = 'sx';
options.COMP      = 1;
% options.TIME      = data.Ts5D;
options.TIME      = 900:10:2000;
data.a = data.EPS * 2000;
create_film(data,options,'.gif')
end

if 0
%% 2D snapshots
% Options
options.INTERP    = 1;
options.POLARPLOT = 0;
options.AXISEQUAL = 1;
% options.NAME      = '\phi';
options.NAME      = 'v_y';
% options.NAME      = 'n_e^{NZ}';
% options.NAME      = '\Gamma_x';
% options.NAME      = 'k^2n_e';
options.PLAN      = 'xy';
% options.NAME      = 'f_e';
% options.PLAN      = 'sx';
options.COMP      = 'avg';
options.TIME      = [000:200:1000];
data.a = data.EPS * 2000;
fig = photomaton(data,options);
save_figure(data,fig)
end

if 0
%% 3D plot on the geometry
options.INTERP    = 1;
options.NAME      = 'n_i';
options.PLANES    = 1;
options.TIME      = [0 500 1000];
data.rho_o_R      = 1e-3; % Sound larmor radius over Machine size ratio
fig = show_geometry(data,options);
save_figure(data,fig)
end

if 0
%% Kinetic distribution function sqrt(<f_a^2>xy) (GENE vsp)
options.SPAR      = linspace(-3,3,64)+(6/127/2);
options.XPERP     = linspace( 0,6,64);
% options.SPAR      = vp';
% options.XPERP     = mu';
options.Z         = 1;
options.T         = 5500;
options.CTR       = 1;
options.ONED      = 0;
fig = plot_fa(data,options);
save_figure(data,fig)
end

if 0
%% Hermite-Laguerre spectrum
% options.TIME = 'avg';
options.TIME = 1000:4000;
options.P2J  = 1;
options.ST   = 0;
options.NORMALIZED = 0;
fig = show_moments_spectrum(data,options);
save_figure(data,fig)
end

if 0
%% 1D spectrum
options.TIME   = 5000:10:5050;
options.NORM   = 1;
options.NAME   = '\phi';
% options.NAME      = 'n_i';
% options.NAME      ='\Gamma_x';
options.PLAN   = 'kxky';
options.COMPZ  = 1;
options.OK     = 0;
options.COMPXY = 'sum';
options.COMPT  = 'avg';
options.PLOT   = 'semilogy';
fig = spectrum_1D(data,options);
% save_figure(data,fig)
end

if 0
%% 1D real plot
options.TIME   = 1500:100:2500;
options.NORM   = 0;
% options.NAME   = '\phi';
% options.NAME      = 'n_i';
% options.NAME      ='\Gamma_x';
options.NAME      ='s_y';
options.COMPZ  = 1;
options.COMPXY = 'avg';
options.COMPT  = 'avg';
options.MOVMT  = 1;
fig = real_plot_1D(data,options);
% save_figure(data,fig)
end

if 0
%% Mode evolution
options.NORMALIZED = 0;
options.K2PLOT = 0.01:0.01:1.0;
options.TIME   = 4900:1:5100;
options.NMA    = 1;
options.NMODES = 20;
fig = mode_growth_meter(data,options);
save_figure(data,fig)
end

if 0
%% ZF caracteristics (space time diagrams
TAVG_0 = 200; TAVG_1 = 3000; % Averaging times duration
% chose your field to plot in spacetime diag (uzf,szf,Gx)
fig = ZF_spacetime(data,TAVG_0,TAVG_1);
save_figure(data,fig)
end
