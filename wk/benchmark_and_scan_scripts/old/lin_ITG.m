%% QUICK RUN SCRIPT
% This script create a directory in /results and run a simulation directly
% from matlab framework. It is meant to run only small problems in linear
% for benchmark and debugging purpose since it makes matlab "busy"
%

gyacomodir  = pwd;
gyacomodir = gyacomodir(1:end-2);
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add% EXECNAME = 'gyacomo_1.0';
% SIMID   = 'linear_CBC_benchmark_GX';  % Name of the simulation
SIMID   = 'dbg';  % Name of the simulation
RUN     = 1; % To run or just to load
default_plots_options
% EXECNAME = 'gyacomo_debug';
EXECNAME = 'gyacomo';
% EXECNAME = 'gyacomo_alpha';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
CLUSTER.TIME  = '99:00:00'; % allocation time hh:mm:ss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 0.0000;           % Collision frequency
TAU     = 1.0;            % e/i temperature ratio
K_Ne    = 0*2.22;            % ele Density '''
K_Te    = 0*6.96;            % ele Temperature '''
K_Ni    = 1*2.22;            % ion Density gradient drive
K_Ti    = 6.96;            % ion Temperature '''
% SIGMA_E = 0.05196152422706632;   % mass ratio sqrt(m_a/m_i) (correct = 0.0233380)
SIGMA_E = 0.0233380;   % mass ratio sqrt(m_a/m_i) (correct = 0.0233380)
KIN_E   = 0;         % 1: kinetic electrons, 2: adiabatic electrons
BETA    = 1e-4;     % electron plasma beta
%% GRID PARAMETERS
P = 6;
J = 2;%P/2;
PMAXE   = P;     % Hermite basis size of electrons
JMAXE   = J;     % Laguerre "
PMAXI   = P;     % " ions
JMAXI   = J;     % "
NX      = 8;     % real space x-gridpoints
NY      = 12;     %     ''     y-gridpoints
LX      = 2*pi/0.1;   % Size of the squared frequency domain
LY      = 2*pi/0.1;   % Size of the squared frequency domain
NZ      = 24;    % number of perpendicular planes (parallel grid)
NPOL    = 1;
SG      = 0;     % Staggered z grids option
NEXC    = 1;     % To extend Lx if needed (Lx = Nexc/(kymin*shear))
%% GEOMETRY
GEOMETRY= 's-alpha';
% GEOMETRY= 'miller';
EPS     = 0.18;   % inverse aspect ratio
Q0      = 1.4;    % safety factor
SHEAR   = 0.8;    % magnetic shear
KAPPA   = 1.0;    % elongation
DELTA   = 0.0;    % triangularity
ZETA    = 0.0;    % squareness
% PARALLEL_BC = 'dirichlet'; %'dirichlet','periodic','shearless','disconnected'
PARALLEL_BC = 'cyclic'; %'dirichlet','periodic','shearless','disconnected'
SHIFT_Y = 0.0;
%% TIME PARMETERS
TMAX    = 25;  % Maximal time unit
% DT      = 1e-2;   % Time step
DT      = 2e-2;   % Time step
SPS0D   = 1;      % Sampling per time unit for 2D arrays
SPS2D   = -1;      % Sampling per time unit for 2D arrays
SPS3D   = 1;      % Sampling per time unit for 2D arrays
SPS5D   = 1/20;    % Sampling per time unit for 5D arrays
SPSCP   = 0;    % Sampling per time unit for checkpoints
JOB2LOAD= -1;
%% OPTIONS
LINEARITY = 'linear';   % activate non-linearity (is cancelled if KXEQ0 = 1)
% Collision operator
% (LB:L.Bernstein, DG:Dougherty, SG:Sugama, LR: Lorentz, LD: Landau)
CO      = 'DG';
GKCO    = 0; % gyrokinetic operator
ABCO    = 1; % INTERSPECIES collisions
INIT_ZF = 0; ZF_AMP = 0.0;
CLOS    = 0;   % Closure model (0: =0 truncation, 1: v^Nmax closure (p+2j<=Pmax))s
NL_CLOS = 0;   % nonlinear closure model (-2:nmax=jmax; -1:nmax=jmax-j; >=0:nmax=NL_CLOS)
KERN    = 0;   % Kernel model (0 : GK)
INIT_OPT= 'phi';   % Start simulation with a noisy mom00/phi/allmom
NUMERICAL_SCHEME = 'RK4'; % RK2,SSPx_RK2,RK3,SSP_RK3,SSPx_RK3,IMEX_SSP2,ARK2,RK4,DOPRI5
%% OUTPUTS
W_DOUBLE = 1;
W_GAMMA  = 1; W_HF     = 1;
W_PHI    = 1; W_NA00   = 1;
W_DENS   = 1; W_TEMP   = 1;
W_NAPJ   = 1; W_SAPJ   = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unused
HD_CO   = 0.0;    % Hyper diffusivity cutoff ratio
MU      = 0.0; % Hyperdiffusivity coefficient
INIT_BLOB = 0; WIPE_TURB = 0; ACT_ON_MODES = 0;
MU_X    = MU;     %
MU_Y    = MU;     %
N_HD    = 4;
MU_Z    = 1.0;     %
HYP_V   = 'dvpar4';
MU_P    = 0.5;     %
MU_J    = 0.0;     %
LAMBDAD = 0.0;
NOISE0  = 1.0e-5; % Init noise amplitude
BCKGD0  = 0.0;    % Init background
k_gB   = 1.0;
k_cB   = 1.0;
COLL_KCUT = 1000;
%%-------------------------------------------------------------------------
%% RUN
setup
% system(['rm fort*.90']);
% Run linear simulation
if RUN
    system(['cd ../results/',SIMID,'/',PARAMS,'/; time mpirun -np 4 ',gyacomodir,'bin/',EXECNAME,' 1 4 1 0; cd ../../../wk'])
%     system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 4 ',gyacomodir,'bin/',EXECNAME,' 1 2 2 0; cd ../../../wk'])
%     system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 2 ',gyacomodir,'bin/',EXECNAME,' 1 2 1 0; cd ../../../wk'])
%     system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 6 ',gyacomodir,'bin/',EXECNAME,' 1 2 3 0; cd ../../../wk'])
%     system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 6 ',gyacomodir,'bin/',EXECNAME,' 1 6 1 0; cd ../../../wk'])
%     system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 1 ',gyacomodir,'bin/',EXECNAME,' 1 1 1 0; cd ../../../wk'])
end

%% Load results
%%
filename = [SIMID,'/',PARAMS,'/'];
LOCALDIR = [gyacomodir,'results/',filename,'/'];
FIGDIR   = LOCALDIR;
% Load outputs from jobnummin up to jobnummax
JOBNUMMIN = 00; JOBNUMMAX = 01;
data = compile_results(LOCALDIR,JOBNUMMIN,JOBNUMMAX); %Compile the results from first output found to JOBNUMMAX if existing
data.FIGDIR = ['../results/',SIMID,'/',PARAMS,'/'];

%% Short analysis
if 0
%% linear growth rate (adapted for 2D zpinch and fluxtube)
options.TRANGE = [0.5 1]*data.Ts3D(end);
options.NPLOTS = 3; % 1 for only growth rate and error, 2 for omega local evolution, 3 for plot according to z
options.GOK    = 0; %plot 0: gamma 1: gamma/k 2: gamma^2/k^3
lg = compute_fluxtube_growth_rate(data,options);
[gmax,     kmax] = max(lg.g_ky(:,end));
[gmaxok, kmaxok] = max(lg.g_ky(:,end)./lg.ky);
msg = sprintf('gmax = %2.2f, kmax = %2.2f',gmax,lg.ky(kmax)); disp(msg);
msg = sprintf('gmax/k = %2.2f, kmax/k = %2.2f',gmaxok,lg.ky(kmaxok)); disp(msg);
end

if 0
%% Ballooning plot
options.time_2_plot = [10 50];
options.kymodes     = [0.1 0.2 0.4];
options.normalized  = 1;
% options.field       = 'phi';
fig = plot_ballooning(data,options);
end


if 0
%% linear growth rate for 3D Zpinch (kz fourier transform)
trange = [0.5 1]*data.Ts3D(end);
options.keq0 = 1; % chose to plot planes at k=0 or max
options.kxky = 1;
options.kzkx = 0;
options.kzky = 0;
[lg, fig] = compute_3D_zpinch_growth_rate(data,trange,options);
save_figure(data,fig)
end
if 1
%% Mode evolution
options.NORMALIZED = 1;
options.TIME   = [000:9000];
options.KX_TW  = [0.5 1]*data.Ts3D(end); %kx Growth rate time window
options.KY_TW  = [0.5 1]*data.Ts3D(end);  %ky Growth rate time window
options.NMA    = 1;
options.NMODES = 800;
options.iz     = 'avg'; % avg or index
options.ik     = 1; % sum, max or index
options.fftz.flag = 0;
fig = mode_growth_meter(data,options);
end


if 0
%% RH TEST
ikx = 2; iky = 2; t0 = 0; t1 = data.Ts3D(end);
[~, it0] = min(abs(t0-data.Ts3D));[~, it1] = min(abs(t1-data.Ts3D));
plt = @(x) squeeze(mean(real(x(iky,ikx,:,it0:it1)),3))./squeeze(mean(real(x(iky,ikx,:,it0)),3));
figure
plot(data.Ts3D(it0:it1), plt(data.PHI));
xlabel('$t$'); ylabel('$\phi_z(t)/\phi_z(0)$')
title(sprintf('$k_x=$%2.2f, $k_y=$%2.2f',data.kx(ikx),data.ky(iky)))
end
