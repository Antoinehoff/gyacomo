addpath(genpath('../matlab')) % ... add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
CLUSTER.TIME  = '99:00:00'; % allocation time hh:mm:ss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 0.0;   % Collision frequency
TAU     = 1.0;    % e/i temperature ratio
ETAB    = 0.6;    % Magnetic gradient
ETAN    = 1.0;    % Density gradient
ETAT    = 0.0;    % Temperature gradient
HD_CO   = 0.5;    % Hyper diffusivity cutoff ratio
NU_HYP  = 0.1;
%% GRID PARAMETERS
N       = 60;     % Frequency gridpoints (Nkr = N/2)
L       = 70;     % Size of the squared frequency domain
PMAXE   = 2;     % Highest electron Hermite polynomial degree
JMAXE   = 1;     % Highest ''       Laguerre ''
PMAXI   = 2;     % Highest ion      Hermite polynomial degree
JMAXI   = 1;     % Highest ''       Laguerre ''
%% TIME PARAMETERS
TMAX    = 60;  % Maximal time unit
DT      = 1e-2;   % Time step
SPS0D   = 1;    % Sampling per time unit for profiler
SPS2D   = 1;      % Sampling per time unit for 2D arrays
SPS5D   = 1;    % Sampling per time unit for 5D arrays
SPSCP   = 1;    % Sampling per time unit for checkpoints/10
RESTART = 1;      % To restart from last checkpoint
JOB2LOAD= 1;
%% OPTIONS
% SIMID   = 'local_nu_%0.0e';  % Name of the simulation
% SIMID   = sprintf(SIMID,NU);
% SIMID   = 'Marconi_DGGK_nu_1e+00';  % Name of the simulation
SIMID   = 'debug_cp_load';  % Name of the simulation
CO      = -3;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty, -3 : GK Dougherty)
CLOS    = 0;   % Closure model (0 : =0 truncation, 1 : n+j = min(nmax,n+j), 2: odd/even adapted)
KERN    = 0;   % Kernel model (0 : GK)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% unused
KR0KH   = 0; A0KH = 0; % Background phi mode to drive Ray-Tay inst.
KREQ0   = 0;      % put kr = 0
KPAR    = 0.0;    % Parellel wave vector component
LAMBDAD = 0.0;
NON_LIN = 1 *(1-KREQ0);   % activate non-linearity (is cancelled if KREQ0 = 1)
kmax    = N*pi/L;% Highest fourier mode
MU      = NU_HYP/(HD_CO*kmax)^4 % Hyperdiffusivity coefficient
NOISE0  = 1.0e-5;

%% Setup and file management
setup
system('rm fort.90');