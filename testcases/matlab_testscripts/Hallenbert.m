addpath(genpath('../matlab')) % ... add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
CLUSTER.TIME  = '99:00:00'; % allocation time hh:mm:ss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 0.1;      % Collision frequency
K_N     = 1.5;       % Density gradient drive
K_T     = 0.25*K_N;  % Temperature '''
K_E     = 0.0;       % Electrostat gradient
SIGMA_E = 0.0233380; % mass ratio sqrt(m_a/m_i) (correct = 0.0233380)
MU_X    = 0.0;      % X hyperdiffusivity
MU_Y    = 0.0;      % Y ''
KIN_E   = 1;         % Kinetic (1) or adiabatic (2) electron model
%% GRID PARAMETERS
NX      = 150;     % Spatial radial resolution ( = 2x radial modes)
LX      = 120;    % Radial window size
NY      = 50;     % Spatial azimuthal resolution (= azim modes)
LY      = 120;    % Azimuthal window size
NZ      = 1;     % number of perpendicular planes (parallel grid)
P       = 4;
J       = 2;
%% GEOMETRY PARAMETERS
Q0      = 1.0;       % safety factor
SHEAR   = 0.0;       % magnetic shear
EPS     = 0.0;      % inverse aspect ratio
k_gB   = 1.0;       % Magnetic  gradient
k_cB   = 1.0;       % Magnetic  curvature
SG      = 0;         % Staggered z grids option
%% TIME PARAMETERS
TMAX    = 2000;  % Maximal time unit
DT      = 1e-3;   % Time step
SPS0D   = 1;      % Sampling per time unit for profiler
SPS2D   = 1;      % Sampling per time unit for 2D arrays
SPS3D   = 1;      % Sampling per time unit for 3D arrays
SPS5D   = 1/100;  % Sampling per time unit for 5D arrays
SPSCP   = 0;    % Sampling per time unit for checkpoints/10
JOB2LOAD= -1;
%% OPTIONS AND NAMING
% Collision operator
% (LB:L.Bernstein, DG:Dougherty, SG:Sugama, LR: Lorentz, LD: Landau)
CO      = 'LD';
GKCO    = 1; % gyrokinetic operator
ABCO    = 1; % INTERSPECIES collisions
NL_CLOS = -1;   % nonlinear closure model (-2: nmax = jmax, -1: nmax = jmax-j, >=0 : nmax = NL_CLOS)
SIMID   = 'Hallenbert_nu_1e-01';  % Name of the simulation
% SIMID   = 'debug';  % Name of the simulation
LINEARITY  = 'nonlinear';   % (nonlinear, semilinear, linear)
% INIT options
INIT_PHI  = 1;   % Start simulation with a noisy phi (0= noisy moments 00)
INIT_ZF   = 0; ZF_AMP = 0.0;
INIT_BLOB = 0; WIPE_TURB = 0; ACT_ON_MODES = 'donothing';
%% OUTPUTS
W_DOUBLE = 1;
W_GAMMA  = 1; W_HF     = 1;
W_PHI    = 1; W_NA00   = 1;
W_DENS   = 1; W_TEMP   = 1;
W_NAPJ   = 1; W_SAPJ   = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% unused
PMAXE   = P;    % Highest electron Hermite polynomial degree
JMAXE   = J;     % Highest ''       Laguerre ''
PMAXI   = P;     % Highest ion      Hermite polynomial degree
JMAXI   = J;     % Highest ''       Laguerre ''
KERN    = 0;   % Kernel model (0 : GK)
KX0KH   = 0; A0KH = 0; % Background phi mode to drive Ray-Tay inst.
KPAR    = 0.0;    % Parellel wave vector component
LAMBDAD = 0.0;
NOISE0  = 1.0e-4;
BCKGD0  = 0;    % Init background
TAU     = 1.0;    % e/i temperature ratio
MU_P    = 0.0;     % Hermite  hyperdiffusivity -mu_p*(d/dvpar)^4 f
MU_J    = 0.0;     % Laguerre hyperdiffusivity -mu_j*(d/dvperp)^4 f
%% Setup and file management
setup
system('rm fort*.90');
outfile = [BASIC.RESDIR,'out.txt'];
disp(outfile);
