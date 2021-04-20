%clear all;
addpath(genpath('../matlab')) % ... add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLUSTER PARAMETERS
CLUSTER.TIME  = '24:00:00'; % allocation time hh:mm:ss
NP_P          = 2;          % MPI processes along p  
NP_KR         = 18;         % MPI processes along kr
CLUSTER.PART  = 'normal';   % debug or normal
if(strcmp(CLUSTER.PART,'debug')); CLUSTER.TIME  = '00:30:00'; end;
CLUSTER.MEM   = '12GB';     % Memory
CLUSTER.JNAME = 'gamma_inf';% Job name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 0.1;   % Collision frequency
ETAB    = 0.6;   % Magnetic gradient
NU_HYP  = 1.0;   % Hyperdiffusivity coefficient
%% GRID PARAMETERS
N       = 200;   % Frequency gridpoints (Nkr = N/2)
L       = 120;   % Size of the squared frequency domain
P       = 12;    % Electron and Ion highest Hermite polynomial degree
J       = 06;    % Electron and Ion highest Laguerre polynomial degree
MU_P    = 0;     % Hermite  hyperdiffusivity -mu_p*(d/dvpar)^4 f
MU_J    = 0;     % Laguerre hyperdiffusivity -mu_j*(d/dvperp)^4 f
%% TIME PARAMETERS
TMAX    = 1000;  % Maximal time unit
DT      = 1e-2;  % Time step
SPS0D   = 1;     % Sampling per time unit for profiler
SPS2D   = 1;     % Sampling per time unit for 2D arrays
SPS5D   = 1/40;  % Sampling per time unit for 5D arrays
SPSCP   = 0;     % Sampling per time unit for checkpoints
RESTART = 1;     % To restart from last checkpoint
JOB2LOAD= 2;
%% OPTIONS
SIMID   = ['HeLaZ_v2.5_eta_',num2str(ETAB),'_nu_%0.0e'];  % Name of the simulation
% SIMID   = 'test_marconi_sugama';  % Name of the simulation
SIMID   = sprintf(SIMID,NU);
PREFIX  =[];
% PREFIX  = sprintf('%d_%d_',NP_P, NP_KR);
% (0 : L.Bernstein, 1 : Dougherty, 2: Sugama, 3 : Full Couloumb ; +/- for GK/DK)
CO      = 1;
CLOS    = 0;   % Closure model (0: =0 truncation, 1: semi coll, 2: Copy closure J+1 = J, P+2 = P)
NL_CLOS = 1;   % nonlinear closure model (0: =0 nmax = jmax, 1: nmax = jmax-j, >1 : nmax = NL_CLOS)
KERN    = 0;   % Kernel model (0 : GK)
INIT_PHI= 1;   % Start simulation with a noisy phi and moments
%% OUTPUTS
W_DOUBLE = 1;
W_GAMMA  = 1;
W_PHI    = 1;
W_NA00   = 1;
W_NAPJ   = 1;
W_SAPJ   = 0;
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
% Compute processes distribution
Ntot = NP_P * NP_KR;
Nnodes = ceil(Ntot/36);
Nppn   = Ntot/Nnodes; 
CLUSTER.NODES =  num2str(Nnodes);  % MPI process along p
CLUSTER.NTPN  =  num2str(Nppn); % MPI process along kr
CLUSTER.CPUPT = '1';        % CPU per task
CLUSTER.NTPC  = '1';        % N tasks per core (openmp threads)
%% Run file management scripts
setup
write_sbash_daint
system('rm fort.90 setup_and_run.sh batch_script.sh');
disp('done');
if(mod(NP_P*NP_KR,36)~= 0)
    disp('WARNING : unused cores (ntot cores must be a 36 multiple)');
end