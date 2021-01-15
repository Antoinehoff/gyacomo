addpath(genpath('../matlab')) % ... add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
CLUSTER.TIME  = '24:00:00'; % allocation time hh:mm:ss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 1e-1;   % Collision frequency
TAU     = 1.0;    % e/i temperature ratio
ETAB    = 0.5;    % Magnetic gradient
ETAN    = 1.0;    % Density gradient
ETAT    = 0.0;    % Temperature gradient
HD_CO   = 0.5;    % Hyper diffusivity cutoff ratio
NU_HYP  = 0.1;
%% GRID PARAMETERS
N       = 128;     % Frequency gridpoints (Nkr = N/2)
L       = 66;     % Size of the squared frequency domain
PMAXE   = 4;     % Highest electron Hermite polynomial degree
JMAXE   = 3;     % Highest ''       Laguerre ''
PMAXI   = 4;     % Highest ion      Hermite polynomial degree
JMAXI   = 3;     % Highest ''       Laguerre ''
%% TIME PARAMETERS
TMAX    = 500;  % Maximal time unit
DT      = 3e-2;   % Time step
SPS0D   = 1/DT;    % Sampling per time unit for profiler
SPS2D   = 1;      % Sampling per time unit for 2D arrays
SPS5D   = 1/5;    % Sampling per time unit for 5D arrays
SPSCP   = 1/10;    % Sampling per time unit for checkpoints
RESTART = 1;      % To restart from last checkpoint
JOB2LOAD= 3;
%% OPTIONS
% SIMID   = 'test_pp2=p_closure';  % Name of the simulation
SIMID   = 'test_truncations';  % Name of the simulation
CO      = -3;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty, -3 : GK Dougherty)
CLOS   = 0;   % Truncation method (0 : =0 closure, 1 : n+j = min(nmax,n+j))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% unused
KR0KH   = 0; A0KH = 0; % Background phi mode to drive Ray-Tay inst.
KREQ0   = 0;      % put kr = 0
KPAR    = 0.0;    % Parellel wave vector component
LAMBDAD = 0.0;
NON_LIN = 1 *(1-KREQ0);   % activate non-linearity (is cancelled if KREQ0 = 1)
LOAD_MARCONI = 0;
kmax    = N*pi/L;% Highest fourier mode
MU      = NU_HYP/(HD_CO*kmax)^4 % Hyperdiffusivity coefficient
NOISE0  = 1.0e-5;

%% Setup and file management
setup
system('rm fort.90');