clear all;
addpath(genpath('../matlab')) % ... add
SUBMIT = 1; % To submit the job automatically
% EXECNAME = 'helaz_dbg';
  EXECNAME = 'helaz_2.62';
for ETAB = [0.6 0.7]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLUSTER PARAMETERS
CLUSTER.PART  = 'prod';     % dbg or prod
CLUSTER.TIME  = '24:00:00'; % allocation time hh:mm:ss
if(strcmp(CLUSTER.PART,'dbg')); CLUSTER.TIME  = '00:30:00'; end;
CLUSTER.MEM   = '128GB';     % Memory
CLUSTER.JNAME = 'HeLaZ';% Job name
NP_P          = 2;          % MPI processes along p  
NP_KR         = 24;         % MPI processes along kr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 1e-1;   % Collision frequency
% ETAB    = 0.7;   % Magnetic gradient
NU_HYP  = 1.0;   % Hyperdiffusivity coefficient
NL_CLOS = -1;   % nonlinear closure model (-2: nmax = jmax, -1: nmax = jmax-j, >=0 : nmax = NL_CLOS)
% (0 : L.Bernstein, 1 : Dougherty, 2: Sugama, 3 : Full Couloumb ; +/- for GK/DK)
CO      = -2;
INIT_ZF = 0; ZF_AMP = 0.0;
%% GRID PARAMETERS
N       = 200;    % Frequency gridpoints (Nkr = N/2)
L       = 120;    % Size of the squared frequency domain
P       = 06;     % Electron and Ion highest Hermite polynomial degree
J       = 03;     % Electron and Ion highest Laguerre polynomial degree
MU_P    = 0.0;% Hermite  hyperdiffusivity -mu_p*(d/dvpar)^4 f
MU_J    = 0.0;% Laguerre hyperdiffusivity -mu_j*(d/dvperp)^4 f
%% TIME PARAMETERS
TMAX    = 5000;  % Maximal time unit
DT      = 1e-2;  % Time step
SPS0D   = 1;     % Sampling per time unit for profiler
SPS2D   = 1/4;     % Sampling per time unit for 2D arrays
SPS5D   = 1/100;  % Sampling per time unit for 5D arrays
SPSCP   = 0;     % Sampling per time unit for checkpoints
RESTART = 0;     % To restart from last checkpoint
JOB2LOAD= 0;
%% Naming
% SIMID   = 'test';  % Name of the simulation
SIMID   = ['v2.6_P_',num2str(P),'_J_',num2str(J)];  % Name of the simulation
PREFIX  =[];
% PREFIX  = sprintf('%d_%d_',NP_P, NP_KR);
%% Options
CLOS    = 0;   % Closure model (0: =0 truncation, 1: semi coll, 2: Copy closure J+1 = J, P+2 = P)
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
% kmax    = 2/3*N*pi/L;% Highest fourier mode with AA
HD_CO   = 0.5;    % Hyper diffusivity cutoff ratio
MU      = NU_HYP/(HD_CO*kmax)^4; % Hyperdiffusivity coefficient
NOISE0  = 1.0e-5;
ETAT    = 0.0;    % Temperature gradient
ETAN    = 1.0;    % Density gradient
TAU     = 1.0;    % e/i temperature ratio
% Compute processes distribution
Ntot = NP_P * NP_KR;
Nnodes = ceil(Ntot/48);
Nppn   = Ntot/Nnodes; 
CLUSTER.NODES =  num2str(Nnodes);  % MPI process along p
CLUSTER.NTPN  =  num2str(Nppn); % MPI process along kr
CLUSTER.CPUPT = '1';        % CPU per task
%% Run file management scripts
setup
write_sbash_marconi
system('rm fort.90 setup_and_run.sh batch_script.sh');
if(mod(NP_P*NP_KR,48)~= 0)
    disp('WARNING : unused cores (ntot cores must be a 48 multiple)');
end
if(SUBMIT)
    system('ssh ahoffman@login.marconi.cineca.it sh HeLaZ/wk/setup_and_run.sh');
end
disp('done');
end