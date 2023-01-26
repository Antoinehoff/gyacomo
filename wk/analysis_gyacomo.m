%% UNCOMMENT FOR TUTORIAL
gyacomodir = pwd; gyacomodir = gyacomodir(1:end-2); % get code directory
% resdir = '.'; %Name of the directory where the results are located
% JOBNUMMIN = 00; JOBNUMMAX = 10;
%%
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add

%% Load the results
LOCALDIR  = [gyacomodir,resdir,'/'];
MISCDIR   = ['/misc/gyacomo_outputs/',resdir,'/']; %For long term storage
system(['mkdir -p ',MISCDIR]);
system(['mkdir -p ',LOCALDIR]);
% CMD = ['rsync ', LOCALDIR,'outputs* ',MISCDIR]; disp(CMD);
% system(CMD);
% Load outputs from jobnummin up to jobnummax
data = compile_results(MISCDIR,JOBNUMMIN,JOBNUMMAX); %Compile the results from first output found to JOBNUMMAX if existing
data.localdir = LOCALDIR;
data.FIGDIR   = LOCALDIR;
data.folder   = LOCALDIR;
data.CODENAME = 'GYACOMO';
%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default_plots_options
disp('Plots')
FMT = '.fig';

if 1
%% Space time diagramm (fig 11 Ivanov 2020)
% data.scale = 1;%/(data.Nx*data.Ny)^2;
i_ = 1; 
disp([num2str(data.TJOB_SE(i_)),' ',num2str(data.TJOB_SE(i_+1))])
disp([num2str(data.NU_EVOL(i_)),' ',num2str(data.NU_EVOL(i_+1))])
options.TAVG_0   = data.TJOB_SE(i_)+600;%0.4*data.Ts3D(end);
options.TAVG_1   = data.TJOB_SE(i_+1);%0.9*data.Ts3D(end); % Averaging times duration
options.NCUT     = 4;              % Number of cuts for averaging and error estimation
options.NMVA     = 1;              % Moving average for time traces
% options.ST_FIELD = '\Gamma_x';   % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.ST_FIELD = '\phi';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.INTERP   = 0;
options.RESOLUTION = 256;
fig = plot_radial_transport_and_spacetime(data,options,'GYACOMO');
% save_figure(data,fig,'.png')
end

if 0
%% statistical transport averaging
gavg =[]; gstd = [];
for i_ = 1:2:numel(data.TJOB_SE) 
% i_ = 3; 
disp([num2str(data.TJOB_SE(i_)),' ',num2str(data.TJOB_SE(i_+1))])
disp([num2str(data.NU_EVOL(i_)),' ',num2str(data.NU_EVOL(i_+1))])
options.T = [data.TJOB_SE(i_) data.TJOB_SE(i_+1)];
options.NPLOTS = 0;
[fig, gavg_, gstd_] = statistical_transport_averaging(data,options);
gavg = [gavg gavg_]; gstd = [gstd gstd_];
end
disp(gavg); disp(gstd);
end
if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
options.INTERP    = 1;
options.POLARPLOT = 0;
options.NAME      = '\phi';
% options.NAME      = '\omega_z';
% options.NAME     = 'N_i^{00}';
% options.NAME      = 'v_x';
% options.NAME      = 'n_i^{NZ}';
% options.NAME      = '\Gamma_x';
% options.NAME      = 'n_i';
% options.NAME      = 'n_i-n_e';
options.PLAN      = 'xy';
options.NAME      = 'f_i';
options.PLAN      = 'sx';
options.COMP      = 'avg';
% options.TIME      = data.Ts5D(end-30:end);
% options.TIME      =  data.Ts3D;
options.TIME      = [0:10000];
data.EPS          = 0.1;
data.a = data.EPS * 2000;
options.RESOLUTION = 256;
create_film(data,options,'.gif')
end

if 0
%% fields snapshots
% Options
options.INTERP    = 0;
options.POLARPLOT = 0;
options.AXISEQUAL = 0;
options.NORMALIZE = 0;
options.NAME      = '\phi';
% options.NAME      = '\psi';
% options.NAME      = '\omega_z';
% options.NAME      = 'n_e';
% options.NAME      = 'n_i-n_e';
% options.NAME      = '\phi^{NZ}';
% options.NAME      = 'N_i^{00}';
% options.NAME      = 'N_i^{00}-N_e^{00}';
% options.NAME      = 's_{Ex}';
% options.NAME      = 'Q_x';
% options.NAME      = 'k^2n_e';
options.PLAN      = 'xy';
options.COMP      = 'avg';
options.TIME      = [500];
options.RESOLUTION = 256;

data.a = data.EPS * 2e3;
fig = photomaton(data,options);
% save_figure(data,fig)
end

if 0
%% plot on the geometry
options.INTERP    = 0;
options.NAME      = '\phi';
options.PLANES    = [1];
options.TIME      = [30];
options.PLT_MTOPO = 1;
options.PLT_FTUBE = 0;
data.EPS = 0.4;
data.rho_o_R      = 3e-3; % Sound larmor radius over Machine size ratio
fig = show_geometry(data,options);
save_figure(data,fig,'.png')
end

if 0
%% Kinetic distribution function sqrt(<f_a^2>xy) (GENE vsp)
% options.SPAR      = linspace(-3,3,32)+(6/127/2);
options.SPAR      = linspace(-3,3,32);
options.XPERP     = linspace( 0,sqrt(6),16).^2;
% options.SPAR      = gene_data.vp';
% options.XPERP     = gene_data.mu';
options.iz        = 'avg';
options.T         = [0.5:0.1:1]*data.Ts3D(end);
% options.PLT_FCT   = 'contour';
% options.PLT_FCT   = 'contourf';
options.PLT_FCT   = 'surfvv';
options.ONED      = 0;
options.non_adiab = 0;
options.SPECIES   = 'i';
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
options.TIME       = [200:500];
options.specie     = 'i';
options.compz      = 'avg';
fig = show_moments_spectrum(data,options);
% fig = show_napjz(data,options);
% save_figure(data,fig,'.png');
end

if 0
%% Time averaged spectrum
options.TIME   = [5000 9000];
options.NORM   =1;
% options.NAME   = '\phi';
% options.NAME      = 'N_i^{00}';
options.NAME   ='\Gamma_x';
options.PLAN   = 'kxky';
options.COMPZ  = 'avg';
options.OK     = 0;
options.COMPXY = 'avg';%'2D'; % avg/sum/max/zero/ 2D plot otherwise
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
options.COMPX  = 1;
options.COMPY  = 2;
options.COMPZ  = 1;
options.COMPT  = 1;
options.MOVMT  = 1;
fig = real_plot_1D(data,options);
% save_figure(data,fig,'.png')
end

if 0
%% Mode evolution
options.NORMALIZED = 1;
options.TIME   = [000:9000];
options.KX_TW  = [25 55]; %kx Growth rate time window
options.KY_TW  = [0 20];  %ky Growth rate time window
options.NMA    = 1;
options.NMODES = 800;
options.iz     = 'avg'; % avg or index
options.ik     = 1; % sum, max or index
options.fftz.flag = 0;
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
options.SHOW_FLUXSURF = 0;
options.SHOW_METRICS  = 1;
fig = plot_metric(data,options);
end

if 0
%% linear growth rate (adapted for 2D zpinch and fluxtube)
options.TRANGE = [0.5 1]*data.Ts3D(end);
options.NPLOTS = 2; % 1 for only growth rate and error, 2 for omega local evolution, 3 for plot according to z
options.GOK    = 0; %plot 0: gamma 1: gamma/k 2: gamma^2/k^3
lg = compute_fluxtube_growth_rate(data,options);
[gmax,     kmax] = max(lg.g_ky(:,end));
[gmaxok, kmaxok] = max(lg.g_ky(:,end)./lg.ky);
msg = sprintf('gmax = %2.2f, kmax = %2.2f',gmax,lg.ky(kmax)); disp(msg);
msg = sprintf('gmax/k = %2.2f, kmax/k = %2.2f',gmaxok,lg.ky(kmaxok)); disp(msg);
end

if 0
%% linear growth rate for 3D Zpinch
trange = [100 200];
options.keq0 = 1; % chose to plot planes at k=0 or max
options.kxky = 1;
options.kzkx = 1;
options.kzky = 1;
options.INTERP = 0;
[lg, fig] = compute_3D_zpinch_growth_rate(data,trange,options);
% save_figure(data,fig,'.png')
end
