% clear all;
addpath(genpath('../matlab')) % ... add
default_plots_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 0e-3;   % Collision frequency
TAU     = 1.0;    % e/i temperature ratio
ETAB    = 0.0;    % Magnetic gradient
ETAN    = 0.0;    % Density gradient
ETAT    = 0.0;    % Temperature gradient
MU      = 1e-4;   % Hyper diffusivity coefficient
LAMBDAD = 0.0; 
NOISE0  = 1.0e-4;
%% GRID PARAMETERS
N       = 64;     % Frequency gridpoints (Nkr = N/2)
L       = 40;     % Size of the squared frequency domain
N0      = 16;     % Periods number for the background sinusoidal profile
KR0KH   = N0*2*pi/L; A0KH = 5/N0; % KH inst.
KREQ0   = 0;      % put kr = 0
PMAXE   = 00;     % Highest electron Hermite polynomial degree
JMAXE   = 07;     % Highest ''       Laguerre ''
PMAXI   = 00;     % Highest ion      Hermite polynomial degree
JMAXI   = 07;     % Highest ''       Laguerre ''
KPAR    = 0.0;    % Parellel wave vector component
%% TIME PARAMETERS 
TMAX    = 20;  % Maximal time unit
DT      = 1e-4;   % Time step
SPS     = 1;      % Sampling per time unit
RESTART = 0;      % To restart from last checkpoint
JOB2LOAD= 00;
%% OPTIONS
SIMID   = 'KH_lin_mode';  % Name of the simulation
NON_LIN = 0 *(1-KREQ0);   % activate non-linearity (is cancelled if KREQ0 = 1)
CO      = 0;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
% DK    = 0;  % Drift kinetic model (put every kernel_n to 0 except n=0 to 1)
NO_E    = 0;  % Remove electrons dynamic
WRITE5D = '.true.';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setup