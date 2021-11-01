addpath(genpath('../matlab')) % ... add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
CLUSTER.TIME  = '99:00:00'; % allocation time hh:mm:ss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 0.1;   % Collision frequency
K_N     = 2.22;      % Density gradient drive
K_T     = 6.0;       % Temperature '''
K_E     = 0.00;    % Electrostat gradient
SIGMA_E = 0.05196;   % mass ratio sqrt(m_a/m_i) (correct = 0.0233380)
NU_HYP  = 0.0;
KIN_E   = 1;         % Kinetic (1) or adiabatic (2) electron model
%% GRID PARAMETERS
NX      = 50;     % Spatial radial resolution ( = 2x radial modes)
LX      = 300;    % Radial window size
NY      = 100;     % Spatial azimuthal resolution (= azim modes)
LY      = 300;    % Azimuthal window size
NZ      = 20;     % number of perpendicular planes (parallel grid)
P       = 4;
J       = 2;
%% GEOMETRY PARAMETERS
Q0      = 2.7;       % safety factor
SHEAR   = 0.0;       % magnetic shear
EPS     = 0.18;      % inverse aspect ratio
GRADB   = 1.0;   % Magnetic  gradient
CURVB   = 1.0;   % Magnetic  curvature
%% TIME PARAMETERS
TMAX    = 10;  % Maximal time unit
DT      = 5e-3;   % Time step
SPS0D   = 1;      % Sampling per time unit for profiler
SPS2D   = 1;      % Sampling per time unit for 2D arrays
SPS3D   = 5;      % Sampling per time unit for 3D arrays
SPS5D   = 1/200;  % Sampling per time unit for 5D arrays
SPSCP   = 0;    % Sampling per time unit for checkpoints/10
JOB2LOAD= -1;
%% OPTIONS AND NAMING
% Collision operator
% (0 : L.Bernstein, 1 : Dougherty, 2: Sugama, 3 : Pitch angle ; 4 : Coulomb; +/- for GK/DK)
CO      = 1;
CLOS    = 0;   % Closure model (0: =0 truncation)
NL_CLOS = 0;   % nonlinear closure model (-2: nmax = jmax, -1: nmax = jmax-j, >=0 : nmax = NL_CLOS)
SIMID   = 'fluxtube_salphaB_s0';  % Name of the simulation
% SIMID   = 'simulation_A';  % Name of the simulation
% SIMID   = ['v3.0_P_',num2str(P),'_J_',num2str(J)];  % Name of the simulation
NON_LIN = 0;   % activate non-linearity (is cancelled if KXEQ0 = 1)
% INIT options
INIT_PHI= 1;   % Start simulation with a noisy phi (0= noisy moments 00)
INIT_ZF   = 0; ZF_AMP = 0.0;
INIT_BLOB = 0; WIPE_TURB = 0; WIPE_ZF = 0;
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
kmax    = NX*pi/LX;% Highest fourier mode
HD_CO   = 0.5;    % Hyper diffusivity cutoff ratio
MU      = NU_HYP/(HD_CO*kmax)^4; % Hyperdiffusivity coefficient
NOISE0  = 1.0e-5;
BCKGD0  = 0.0;    % Init background
TAU     = 1.0;    % e/i temperature ratio
MU_P    = 0.0;     % Hermite  hyperdiffusivity -mu_p*(d/dvpar)^4 f
MU_J    = 0.0;     % Laguerre hyperdiffusivity -mu_j*(d/dvperp)^4 f
%% Setup and file management
setup
system('rm fort*.90');
outfile = [BASIC.RESDIR,'out.txt'];
disp(outfile);
