%% QUICK RUN SCRIPT
% This script creates a directory in /results and runs a simulation directly
% from the Matlab framework. It is meant to run only small problems in linear
% for benchmarking and debugging purposes since it makes Matlab "busy".

%% Set up the paths for the necessary Matlab modules
gyacomodir = pwd;
gyacomodir = gyacomodir(1:end-2);
addpath(genpath([gyacomodir,'matlab'])) % Add matlab module
addpath(genpath([gyacomodir,'matlab/plot'])) % Add plot module
addpath(genpath([gyacomodir,'matlab/compute'])) % Add compute module
addpath(genpath([gyacomodir,'matlab/load'])) % Add load module

%% Set simulation parameters
SIMID   = 'ETG_adiab_i';  % Name of the simulation
RUN = 1; % To run or just to load
default_plots_options
EXECNAME = 'gyacomo23_sp'; % single precision
% EXECNAME = 'gyacomo23_dp'; % double precision
% EXECNAME = 'gyacomo23_debug'; % single precision

%% Set up physical parameters
CLUSTER.TIME = '99:00:00';  % Allocation time hh:mm:ss
NU = 0.00;                   % Collision frequency
TAU = 1.0;                  % e/i temperature ratio
K_Ne    = 2.5;             % ele Density '''
K_Te    = K_Ne/4;             % ele Temperature '''
K_Ni    = 2.5;                % ion Density gradient drive
K_Ti    = K_Ni/4;                % ion Temperature '''
SIGMA_E = 0.0233380;        % mass ratio sqrt(m_a/m_i) (correct = 0.0233380)
NA = 2;                     % number of kinetic species
ADIAB_E = 0;          % adiabatic electron model
ADIAB_I = 0;          % adiabatic ion model
BETA    = 0.000;             % electron plasma beta
MHD_PD  = 0;                % MHD pressure drift
%% Set up grid parameters
P = 4;
J = P/2;%P/2;
PMAX = P;                   % Hermite basis size
JMAX = J;                   % Laguerre basis size
NX = 2;                     % real space x-gridpoints
NY = 12;                    % real space y-gridpoints
LX = 2*pi/0.05;              % Size of the squared frequency domain in x direction
LY = 2*pi/0.2;              % Size of the squared frequency domain in y direction
NZ = 16;                    % number of perpendicular planes (parallel grid)
SG = 0;                     % Staggered z grids option
NEXC = 1;                   % To extend Lx if needed (Lx = Nexc/(kymin*shear))
%% GEOMETRY
% GEOMETRY= 's-alpha';
% GEOMETRY= 'miller';
GEOMETRY= 'z-pinch';
EPS     = 0.18;   % inverse aspect ratio
Q0      = 1.4;    % safety factor
SHEAR   = 0.8;    % magnetic shear
KAPPA   = 1.0;   % elongation
S_KAPPA = 0.0;
DELTA   = 0.0;  % triangularity
S_DELTA = 0.0;
ZETA    = 0.0; % squareness
S_ZETA  = 0.0;
PARALLEL_BC = 'dirichlet'; % Boundary condition for parallel direction ('dirichlet','periodic','shearless','disconnected')
SHIFT_Y = 0.0;    % Shift in the periodic BC in z
NPOL    = 1;       % Number of poloidal turns
PB_PHASE= 0;
%% TIME PARAMETERS
TMAX     = 100;  % Maximal time unit
DT       = 1e-3;   % Time step
DTSAVE0D = 1;      % Sampling per time unit for 0D arrays
DTSAVE2D = -1;     % Sampling per time unit for 2D arrays
DTSAVE3D = 1;      % Sampling per time unit for 3D arrays
DTSAVE5D = 100;     % Sampling per time unit for 5D arrays
JOB2LOAD = -1;     % Start a new simulation serie

%% OPTIONS
LINEARITY = 'linear';   % activate non-linearity (is cancelled if KXEQ0 = 1)
CO        = 'DG';       % Collision operator (LB:L.Bernstein, DG:Dougherty, SG:Sugama, LR: Lorentz, LD: Landau)
GKCO      = 0;          % Gyrokinetic operator
ABCO      = 1;          % INTERSPECIES collisions
INIT_ZF   = 0;          % Initialize zero-field quantities
HRCY_CLOS = 'truncation';   % Closure model for higher order moments
DMAX      = -1;
NLIN_CLOS = 'truncation';   % Nonlinear closure model for higher order moments
NMAX      = 0;
KERN      = 0;   % Kernel model (0 : GK)
INIT_OPT  = 'phi';   % Start simulation with a noisy mom00/phi/allmom
NUMERICAL_SCHEME = 'RK4'; % Numerical integration scheme (RK2,SSPx_RK2,RK3,SSP_RK3,SSPx_RK3,IMEX_SSP2,ARK2,RK4,DOPRI5)

%% OUTPUTS
W_DOUBLE = 1;     % Output flag for double moments
W_GAMMA  = 1;     % Output flag for gamma (Gyrokinetic Energy)
W_HF     = 1;     % Output flag for high-frequency potential energy
W_PHI    = 1;     % Output flag for potential
W_NA00   = 1;     % Output flag for nalpha00 (density of species alpha)
W_DENS   = 1;     % Output flag for total density
W_TEMP   = 1;     % Output flag for temperature
W_NAPJ   = 1;     % Output flag for nalphaparallel (parallel momentum of species alpha)
W_SAPJ   = 0;     % Output flag for saparallel (parallel current of species alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UNUSED PARAMETERS
% These parameters are usually not to play with in linear runs
MU      = 0.0;    % Hyperdiffusivity coefficient
MU_X    = MU;     % Hyperdiffusivity coefficient in x direction
MU_Y    = MU;     % Hyperdiffusivity coefficient in y direction
N_HD    = 4;      % Degree of spatial-hyperdiffusivity
MU_Z    = 2.0;    % Hyperdiffusivity coefficient in z direction
HYP_V   = 'hypcoll'; % Kinetic-hyperdiffusivity model
MU_P    = 0.0;    % Hyperdiffusivity coefficient for Hermite
MU_J    = 0.0;    % Hyperdiffusivity coefficient for Laguerre
LAMBDAD = 0.0;    % Lambda Debye
NOISE0  = 1.0e-4; % Initial noise amplitude
BCKGD0  = 0.0e-8;    % Initial background
k_gB   = 1.0;     % Magnetic gradient strength
k_cB   = 1.0;     % Magnetic curvature strength
COLL_KCUT = 1; % Cutoff for collision operator

%%-------------------------------------------------------------------------
%% RUN
setup
% system(['rm fort*.90']);
% Run linear simulation
if RUN
    MVIN =['cd ../results/',SIMID,'/',PARAMS,'/;'];
    % RUN  =['time mpirun -np 2 ',gyacomodir,'bin/',EXECNAME,' 1 2 1 0;'];
    RUN  =['time mpirun -np 4 ',gyacomodir,'bin/',EXECNAME,' 1 4 1 0;'];
    % RUN  =['time mpirun -np 6 ',gyacomodir,'bin/',EXECNAME,' 1 3 2 0;'];
    % RUN  =['time mpirun -np 1 ',gyacomodir,'bin/',EXECNAME,' 1 1 1 0;'];
    MVOUT='cd ../../../wk;';
    system([MVIN,RUN,MVOUT]);
end

%% Analysis
% load
filename = [SIMID,'/',PARAMS,'/']; % Create the filename based on SIMID and PARAMS
LOCALDIR = [gyacomodir,'results/',filename,'/']; % Create the local directory path based on gyacomodir, results directory, and filename
FIGDIR   = LOCALDIR; % Set FIGDIR to the same path as LOCALDIR
% Load outputs from jobnummin up to jobnummax
J0 = 0; J1 = 0;
data = {}; % Initialize data as an empty cell array
% load grids, inputs, and time traces
data = compile_results_low_mem(data,LOCALDIR,J0,J1); 


if 1 % Activate or not
%% plot mode evolution and growth rates
% Load phi
[data.PHI, data.Ts3D] = compile_results_3D(LOCALDIR,J0,J1,'phi');
options.NORMALIZED = 0; 
options.TIME   = data.Ts3D;
 % Time window to measure the growth of kx/ky modes
options.KX_TW  = [0.2 1]*data.Ts3D(end);
options.KY_TW  = [0.2 1]*data.Ts3D(end);
options.NMA    = 1; % Set NMA option to 1
options.NMODES = 999; % Set how much modes we study
options.iz     = 'avg'; % Compressing z
options.ik     = 1; %
options.fftz.flag = 0; % Set fftz.flag option to 0
fig = mode_growth_meter(data,options); % Call the function mode_growth_meter with data and options as input arguments, and store the result in fig
end


