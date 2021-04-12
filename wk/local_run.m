addpath(genpath('../matlab')) % ... add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
CLUSTER.TIME  = '99:00:00'; % allocation time hh:mm:ss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 1.0;   % Collision frequency
ETAB    = 0.5;    % Magnetic gradient
ETAN    = 1.0;    % Density gradient
NU_HYP  = 1.0;
%% GRID PARAMETERS
N       = 200;     % Frequency gridpoints (Nkr = N/2)
L       = 120;     % Size of the squared frequency domain
PMAXE   = 10;     % Highest electron Hermite polynomial degree
JMAXE   = 05;     % Highest ''       Laguerre ''
PMAXI   = 10;     % Highest ion      Hermite polynomial degree
JMAXI   = 05;     % Highest ''       Laguerre ''
%% TIME PARAMETERS
TMAX    = 50;  % Maximal time unit
DT      = 2e-2;   % Time step
SPS0D   = 1;    % Sampling per time unit for profiler
SPS2D   = 1/2;      % Sampling per time unit for 2D arrays
SPS5D   = 1/4;    % Sampling per time unit for 5D arrays
SPSCP   = 0;    % Sampling per time unit for checkpoints/10
RESTART = 0;      % To restart from last checkpoint
JOB2LOAD= 0;
%% OPTIONS AND NAMING
% Collision operator
% (0 : L.Bernstein, 1 : Dougherty, 2: Sugama, 3 : Full Couloumb ; +/- for GK/DK)
CO      = 2 ;
CLOS    = 0;   % Closure model (0: =0 truncation, 1: semi coll, 2: Copy closure J+1 = J, P+2 = P)
NL_CLOS = -1;   % nonlinear closure model (-2: nmax = jmax, -1: nmax = jmax-j, >=0 : nmax = NL_CLOS)
SIMID   = 'test_SGGK_full_size';  % Name of the simulation
NON_LIN = 1 *(1-KREQ0);   % activate non-linearity (is cancelled if KREQ0 = 1)
%% OUTPUTS
W_DOUBLE = 0;
W_GAMMA  = 1;
W_PHI    = 1;
W_NA00   = 1;
W_NAPJ   = 1;
W_SAPJ   = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% unused
KERN    = 0;   % Kernel model (0 : GK)
KR0KH   = 0; A0KH = 0; % Background phi mode to drive Ray-Tay inst.
KREQ0   = 0;      % put kr = 0
KPAR    = 0.0;    % Parellel wave vector component
LAMBDAD = 0.0;
kmax    = N*pi/L;% Highest fourier mode
HD_CO   = 0.5;    % Hyper diffusivity cutoff ratio
% kmaxcut = 2.5;
MU      = NU_HYP/(HD_CO*kmax)^4 % Hyperdiffusivity coefficient
NOISE0  = 1.0e-5;
TAU     = 1.0;    % e/i temperature ratio
ETAT    = 0.0;    % Temperature gradient
MU_P    = 0.0/PMAXI^2;     % Hermite  hyperdiffusivity -mu_p*(d/dvpar)^4 f
MU_J    = 0.0/JMAXI^3;     % Laguerre hyperdiffusivity -mu_j*(d/dvperp)^4 f
INIT_PHI= 1;   % Start simulation with a noisy phi and moments
%% Setup and file management
setup
system('rm fort.90');