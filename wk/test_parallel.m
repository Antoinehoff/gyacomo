clear all;
addpath(genpath('../matlab')) % ... add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 1e-1;   % Collision frequency
TAU     = 1.0;    % e/i temperature ratio
ETAB    = 0.6;    % Magnetic gradient
ETAN    = 1.0;    % Density gradient
ETAT    = 0.0;    % Temperature gradient
MU      = 5e-4;   % Hyper diffusivity coefficient
NOISE0  = 1.0e-5;
%% GRID PARAMETERS
N       = 400;     % Frequency gridpoints (Nkr = N/2)
L       = 100;     % Size of the squared frequency domain
PMAXE   = 2;     % Highest electron Hermite polynomial degree
JMAXE   = 1;     % Highest ''       Laguerre ''
PMAXI   = 2;     % Highest ion      Hermite polynomial degree
JMAXI   = 1;     % Highest ''       Laguerre ''
%% TIME PARAMETERS 
TMAX    = 200;  % Maximal time unit
DT      = 5e-3;   % Time step
SPS0D   = 1;    % Sampling per time unit for profiler
SPS2D   = 1;      % Sampling per time unit for 2D arrays
SPS5D   = 0.1;    % Sampling per time unit for 5D arrays
RESTART = 1;      % To restart from last checkpoint
JOB2LOAD= 0;
%% OPTIONS
SIMID   = 'test_parallel';  % Name of the simulation
CO      = -1;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% unused
KR0KH   = 0; A0KH = 0; % Background phi mode to drive Ray-Tay inst.
NO_E    = 0;  % Remove electrons dynamic
% DK    = 0;  % Drift kinetic model (put every kernel_n to 0 except n=0 to 1)
KREQ0   = 0;      % put kr = 0
KPAR    = 0.0;    % Parellel wave vector component
LAMBDAD = 0.0; 
NON_LIN = 1 *(1-KREQ0);   % activate non-linearity (is cancelled if KREQ0 = 1)
CANCEL_ODD_P = 0;% Cancels the odd polynomials degree

setup
