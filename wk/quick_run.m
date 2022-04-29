%% QUICK RUN SCRIPT
% This script create a directory in /results and run a simulation directly
% from matlab framework. It is meant to run only small problems in linear
% for benchmark and debugging purpose since it makes matlab "busy"
%
RUN = 1; % To run or just to load
addpath(genpath('../matlab')) % ... add
default_plots_options
HELAZDIR = '/home/ahoffman/HeLaZ/';
EXECNAME = 'helaz3';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
CLUSTER.TIME  = '99:00:00'; % allocation time hh:mm:ss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 0.1;   % Collision frequency
TAU     = 1.0;    % e/i temperature ratio
K_N     = 20;%1.9;%2.22;   % Density gradient drive
K_T     = 20;%0.25*K_N;   % Temperature '''
K_E     = 0.0;   % Electrostat '''
SIGMA_E = 0.05;%0.0233380;   % mass ratio sqrt(m_a/m_i) (correct = 0.0233380)
KIN_E   = 1;     % 1: kinetic electrons, 2: adiabatic electrons
%% GRID PARAMETERS
PMAXE   = 2;     % Hermite basis size of electrons
JMAXE   = 1;     % Laguerre "
PMAXI   = 2;     % " ions
JMAXI   = 1;     % "
NX      = 32;    % real space x-gridpoints
NY      = 32;     %     ''     y-gridpoints
LX      = 300;   % Size of the squared frequency domain
LY      = 300;     % Size of the squared frequency domain
NZ      = 16;     % number of perpendicular planes (parallel grid)
NPOL    = 1;
SG      = 0;     % Staggered z grids option
%% GEOMETRY
% GEOMETRY= 'Z-pinch'; % Z-pinch overwrites q0, shear and eps
GEOMETRY= 's-alpha';
Q0      = 2.5;    % safety factor
SHEAR   = 0.0;    % magnetic shear (Not implemented yet)
EPS     = 0.18;    % inverse aspect ratio
%% TIME PARMETERS
TMAX    = 50;  % Maximal time unit
DT      = 1e-3;   % Time step
SPS0D   = 1;      % Sampling per time unit for 2D arrays
SPS2D   = 0;      % Sampling per time unit for 2D arrays
SPS3D   = 1;      % Sampling per time unit for 2D arrays
SPS5D   = 1;    % Sampling per time unit for 5D arrays
SPSCP   = 0;    % Sampling per time unit for checkpoints
JOB2LOAD= -1;
%% OPTIONS
SIMID   = 'quick_run';  % Name of the simulation
LINEARITY = 'linear';   % activate non-linearity (is cancelled if KXEQ0 = 1)
% Collision operator
% (LB:L.Bernstein, DG:Dougherty, SG:Sugama, LR: Lorentz, LD: Landau)
CO      = 'DG';
GKCO    = 1; % gyrokinetic operator
ABCO    = 1; % interspecies collisions
INIT_ZF = 0; ZF_AMP = 0.0;
CLOS    = 0;   % Closure model (0: =0 truncation, 1: gyrofluid closure (p+2j<=Pmax))s
NL_CLOS = 0;   % nonlinear closure model (-2:nmax=jmax; -1:nmax=jmax-j; >=0:nmax=NL_CLOS)
KERN    = 0;   % Kernel model (0 : GK)
INIT_OPT= 'phi';   % Start simulation with a noisy mom00/phi/allmom
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
kmax    = NX*pi/LX;% Highest fourier mode
MU      = 0.0; % Hyperdiffusivity coefficient
INIT_BLOB = 0; WIPE_TURB = 0; ACT_ON_MODES = 0;
MU_X    = MU;     % 
MU_Y    = MU;     % 
MU_Z    = 0.0;     %
MU_P    = 0.0;     %
MU_J    = 0.0;     % 
LAMBDAD = 0.0;
NOISE0  = 1.0e-5; % Init noise amplitude
BCKGD0  = 0.0;    % Init background
GRADB   = 1.0;
CURVB   = 1.0;
%%-------------------------------------------------------------------------
%% RUN
setup
system(['rm fort*.90']);
% Run linear simulation
if RUN
%     system(['cd ../results/',SIMID,'/',PARAMS,'/; time mpirun -np 1 ',HELAZDIR,'bin/',EXECNAME,' 1 1 1 0; cd ../../../wk'])
%     system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 2 ',HELAZDIR,'bin/',EXECNAME,' 2 1 1 0; cd ../../../wk'])
%     system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 4 ',HELAZDIR,'bin/',EXECNAME,' 1 2 2 0; cd ../../../wk'])
    system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 6 ',HELAZDIR,'bin/',EXECNAME,' 1 3 2 0; cd ../../../wk'])
end

%% Load results
%%
filename = [SIMID,'/',PARAMS,'/'];
LOCALDIR  = [HELAZDIR,'results/',filename,'/'];
% Load outputs from jobnummin up to jobnummax
JOBNUMMIN = 00; JOBNUMMAX = 00; 
data = compile_results(LOCALDIR,JOBNUMMIN,JOBNUMMAX); %Compile the results from first output found to JOBNUMMAX if existing

%% Short analysis
if 1
%% linear growth rate (adapted for 2D zpinch and fluxtube)
trange = [0.5 1]*data.Ts3D(end);
nplots = 1;
lg = compute_fluxtube_growth_rate(data,trange,nplots);
end

if 0
%% Ballooning plot
options.time_2_plot = data.Ts3D(end);
options.normalized  = 1;
options.field       = 'phi';
fig = plot_ballooning(data,options);
end

if 1
%% linear growth rate for 3D Zpinch (kz fourier transform)
trange = [0.5 1]*data.Ts3D(end);
options.keq0 = 1; % chose to plot planes at k=0 or max
options.kxky = 1;
options.kzkx = 0;
options.kzky = 0;
[lg, fig] = compute_3D_zpinch_growth_rate(data,trange,options);
save_figure(data,fig)
end
