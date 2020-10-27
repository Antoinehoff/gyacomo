clear all;
addpath(genpath('../matlab')) % ... add
default_plots_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 1e-3;   % Collision frequency
TAU     = 1.0;    % e/i temperature ratio
ETAB    = 0.5;    % Magnetic gradient
ETAN    = 1.0;    % Density gradient
ETAT    = 0.0;    % Temperature gradient
MU      = 1e-4;   % Hyper diffusivity coefficient
LAMBDAD = 0.0; 
NOISE0  = 5.0e-5;
%% GRID PARAMETERS
N       = 256;     % Frequency gridpoints (Nkr = N/2)
L       = 40;     % Size of the squared frequency domain
KREQ0   = 0;      % put kr = 0
PMAXE   = 02;     % Highest electron Hermite polynomial degree
JMAXE   = 01;     % Highest ''       Laguerre ''
PMAXI   = 02;     % Highest ion      Hermite polynomial degree
JMAXI   = 01;     % Highest ''       Laguerre ''
KPAR    = 0.0;    % Parellel wave vector component
%% TIME PARAMETERS 
TMAX    = 150;  % Maximal time unit
DT    = 1e-2;   % Time step
SPS2D   = 5;      % Sampling per time unit for 2D arrays
SPS5D   = 0.5;    % Sampling per time unit for 5D arrays
RESTART = 0;      % To restart from last checkpoint
JOB2LOAD= 00;
%% OPTIONS
SIMID   = 'ZP_forced_sym';  % Name of the simulation
NON_LIN = 1 *(1-KREQ0);   % activate non-linearity (is cancelled if KREQ0 = 1)
CO      = 0;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unused
KR0KH   = 0; A0KH = 0; % Background phi mode to drive Ray-Tay inst.
NO_E    = 0;  % Remove electrons dynamic
% DK    = 0;  % Drift kinetic model (put every kernel_n to 0 except n=0 to 1)


setup
