addpath(genpath('../matlab')) % ... add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
CLUSTER.TIME  = '99:00:00'; % allocation time hh:mm:ss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 0.1;   % Collision frequency
ETAN    = 1.0/0.6;    % Density gradient drive (R/Ln)
NU_HYP  = 1.0;
%% GRID AND GEOMETRY PARAMETERS
N       = 50;     % Frequency gridpoints (Nkx = N/2)
L       = 30;     % Size of the squared frequency domain
Nz      = 10;      % number of perpendicular planes (parallel grid)
q0      = 1.0;    % safety factor
shear   = 0.0;    % magnetic shear
eps     = 0.0;    % inverse aspect ratio
P       = 2;
J       = 1;
%% TIME PARAMETERS
TMAX    = 500;  % Maximal time unit
DT      = 1e-2;   % Time step
SPS0D   = 1;      % Sampling per time unit for profiler
SPS2D   = 1;      % Sampling per time unit for 2D arrays
SPS3D   = 2;      % Sampling per time unit for 3D arrays
SPS5D   = 1;  % Sampling per time unit for 5D arrays
SPSCP   = 0;    % Sampling per time unit for checkpoints/10
RESTART = 0;      % To restart from last checkpoint
JOB2LOAD= 0;
%% OPTIONS AND NAMING
% Collision operator
% (0 : L.Bernstein, 1 : Dougherty, 2: Sugama, 3 : Pitch angle ; +/- for GK/DK)
CO      = 1;
CLOS    = 0;   % Closure model (0: =0 truncation)
NL_CLOS = 0;   % nonlinear closure model (-2: nmax = jmax, -1: nmax = jmax-j, >=0 : nmax = NL_CLOS)
% SIMID   = 'HD_study';  % Name of the simulation
SIMID   = 'test_3D';  % Name of the simulation
% SIMID   = ['v3.0_P_',num2str(P),'_J_',num2str(J)];  % Name of the simulation
NON_LIN = 1;   % activate non-linearity (is cancelled if KXEQ0 = 1)
% INIT options
INIT_ZF = 0; ZF_AMP = 0.0;
INIT_BLOB = 0; WIPE_TURB = 0;
%% OUTPUTS
W_DOUBLE = 0;
W_GAMMA  = 1;
W_PHI    = 1;
W_NA00   = 1;
W_NAPJ   = 1;
W_SAPJ   = 0;
W_DENS   = 1;
W_TEMP   = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% unused
PMAXE   = P;    % Highest electron Hermite polynomial degree
JMAXE   = J;     % Highest ''       Laguerre ''
PMAXI   = P;     % Highest ion      Hermite polynomial degree
JMAXI   = J;     % Highest ''       Laguerre ''
KERN    = 0;   % Kernel model (0 : GK)
KX0KH   = 0; A0KH = 0; % Background phi mode to drive Ray-Tay inst.
KXEQ0   = 0;      % put kx = 0
KPAR    = 0.0;    % Parellel wave vector component
LAMBDAD = 0.0;
kmax    = N*pi/L;% Highest fourier mode
HD_CO   = 0.5;    % Hyper diffusivity cutoff ratio
% kmaxcut = 2.5;
MU      = NU_HYP/(HD_CO*kmax)^4 % Hyperdiffusivity coefficient
NOISE0  = 1.0e-5;
TAU     = 1.0;    % e/i temperature ratio
ETAT    = 0.0;    % Temperature gradient
ETAB    = 1.0;    % Magnetic gradient (1.0 to set R=LB)
INIT_PHI= 1;   % Start simulation with a noisy phi and moments
MU_P    = 0.0;     % Hermite  hyperdiffusivity -mu_p*(d/dvpar)^4 f
MU_J    = 0.0;     % Laguerre hyperdiffusivity -mu_j*(d/dvperp)^4 f
%% Setup and file management
setup
system('rm fort.90');
outfile = [BASIC.RESDIR,'out.txt'];
disp(outfile);
