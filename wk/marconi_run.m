%clear all;
addpath(genpath('../matlab')) % ... add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLUSTER PARAMETERS
CLUSTER.TIME  = '24:00:00'; % allocation time hh:mm:ss
CLUSTER.NODES = '1';        % MPI process
CLUSTER.CPUPT = '1';        % CPU per task
CLUSTER.NTPN  = '24';       % N tasks per node
CLUSTER.PART  = 'prod';     % dbg or prod
CLUSTER.MEM   = '16GB';     % Memory
CLUSTER.JNAME = 'gamma_inf';     % Job name
%% PHYSICAL PARAMETERS
NU      = 1e-1;   % Collision frequency
TAU     = 1.0;    % e/i temperature ratio
ETAB    = 0.66;    % Magnetic gradient
ETAN    = 1.0;    % Density gradient
ETAT    = 0.0;    % Temperature gradient
HD_CO   = 0.5;    % Hyper diffusivity cutoff ratio
NOISE0  = 1.0e-5;
%% GRID PARAMETERS
N       = 200;     % Frequency gridpoints (Nkr = N/2)
L       = 100;     % Size of the squared frequency domain
P       = 4;       % Electron and Ion highest Hermite polynomial degree
J       = 2;       % Electron and Ion highest Laguerre polynomial degree
%% TIME PARAMETERS
TMAX    = 500;  % Maximal time unit
DT      = 1e-2;   % Time step
SPS0D   = 10;    % Sampling per time unit for profiler
SPS2D   = 1;      % Sampling per time unit for 2D arrays
SPS5D   = 1/10;    % Sampling per time unit for 5D arrays
SPSCP   = 1/10;    % Sampling per time unit for checkpoints
RESTART = 0;      % To restart from last checkpoint
JOB2LOAD= 0;
%% OPTIONS
SIMID   = 'Marconi';  % Name of the simulation
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
PMAXE   = P;    % Highest electron Hermite polynomial degree
JMAXE   = J;     % Highest ''       Laguerre ''
PMAXI   = P;     % Highest ion      Hermite polynomial degree
JMAXI   = J;     % Highest ''       Laguerre ''
kmax    = N*pi/L;% Highest fourier mode
MU      = 0.1/(HD_CO*kmax)^4 % Hyperdiffusivity coefficient
%% Run following scripts
setup

write_sbash

system(['scp {fort.90,batch_script.sh,setup_and_run.sh}',...
    ' ahoffman@login.marconi.cineca.it:/marconi/home/userexternal/ahoffman/HeLaZ/wk']);
LOAD_MARCONI = 1;
disp('done');
