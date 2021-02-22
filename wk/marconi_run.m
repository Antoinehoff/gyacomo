%clear all;
addpath(genpath('../matlab')) % ... add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLUSTER PARAMETERS
CLUSTER.TIME  = '12:00:00'; % allocation time hh:mm:ss
CLUSTER.NODES = '1';        % MPI process
CLUSTER.CPUPT = '1';        % CPU per task
CLUSTER.NTPN  = '24';       % N tasks per node
CLUSTER.PART  = 'prod';     % dbg or prod
CLUSTER.MEM   = '16GB';     % Memory
CLUSTER.JNAME = 'gamma_inf';% Job name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 1.0;   % Collision frequency
ETAB    = 0.6;   % Magnetic gradient
NU_HYP  = 0.1;   % Hyperdiffusivity coefficient
%% GRID PARAMETERS
N       = 150;     % Frequency gridpoints (Nkr = N/2)
L       = 70;     % Size of the squared frequency domain
P       = 10;       % Electron and Ion highest Hermite polynomial degree
J       = 5;       % Electron and Ion highest Laguerre polynomial degree
MU_P    = 0;     % Hermite  hyperdiffusivity -mu_p*(d/dvpar)^4 f
MU_J    = 0;     % Laguerre hyperdiffusivity -mu_j*(d/dvperp)^4 f
%% TIME PARAMETERS
TMAX    = 500;  % Maximal time unit
DT      = 1e-2;  % Time step
SPS0D   = 1;      % Sampling per time unit for profiler
SPS2D   = 1/2;   % Sampling per time unit for 2D arrays
SPS5D   = 1/10;  % Sampling per time unit for 5D arrays
SPSCP   = 1/10;  % Sampling per time unit for checkpoints
RESTART = 1;     % To restart from last checkpoint
JOB2LOAD= 0;
%% OPTIONS
% SIMID   = 'Marconi_DGGK_nu_%0.0e';  % Name of the simulation
SIMID   = ['Marconi_DGGK_eta_',num2str(ETAB),'_nu_%0.0e'];  % Name of the simulation
% SIMID   = 'Marconi_restart';  % Name of the simulation
SIMID   = sprintf(SIMID,NU);
CO      = -3;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty, -3 : GK Dougherty)
CLOS    = 0;   % Closure model (0: =0 truncation, 1: semi coll, 2: Copy closure J+1 = J, P+2 = P)
KERN    = 0;   % Kernel model (0 : GK)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fixed parameters (for current study)
KR0KH   = 0; A0KH = 0; % Background phi mode to drive Ray-Tay inst.
KREQ0   = 0;      % put kr = 0
KPAR    = 0.0;    % Parellel wave vector component
LAMBDAD = 0.0;
NON_LIN = 1 *(1-KREQ0);   % activate non-linearity (is cancelled if KREQ0 = 1)
PMAXE   = P;    % Highest electron Hermite polynomial degree
JMAXE   = J;     % Highest ''       Laguerre ''
PMAXI   = P;     % Highest ion      Hermite polynomial degree
JMAXI   = J;     % Highest ''       Laguerre ''
kmax    = N*pi/L;% Highest fourier mode
HD_CO   = 0.5;    % Hyper diffusivity cutoff ratio
MU      = NU_HYP/(HD_CO*kmax)^4 % Hyperdiffusivity coefficient
NOISE0  = 1.0e-5;
ETAT    = 0.0;    % Temperature gradient
ETAN    = 1.0;    % Density gradient
TAU     = 1.0;    % e/i temperature ratio
%% Run file management scripts
setup
write_sbash
system('rm fort.90 setup_and_run.sh batch_script.sh');
disp('done');
