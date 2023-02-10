clear all;
addpath(genpath('../matlab')) % ... add
SUBMIT = 1; % To submit the job automatically
CHAIN  = 0; % To chain jobs (CHAIN = n will launch n jobs in chain)
% EXECNAME = 'helaz3_dbg';
  EXECNAME = 'helaz3';
  SIMID = 'HP_fig2b_conv';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLUSTER PARAMETERS
CLUSTER.PART  = 'prod';     % dbg or prod
% CLUSTER.PART  = 'dbg';
CLUSTER.TIME  = '24:00:00'; % allocation time hh:mm:ss
if(strcmp(CLUSTER.PART,'dbg')); CLUSTER.TIME  = '00:30:00'; end;
CLUSTER.MEM   = '16GB';     % Memory
CLUSTER.JNAME = 'nu0_b_conv';% Job name
NP_P          = 1;          % MPI processes along p
NP_KX         = 12;         % MPI processes along kx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 0.0;   % Collision frequency
K_N     = 2.5;    % Density gradient drive (R/Ln)
K_T     = K_N/4;    % Temperature gradient
MU      = 0.1;
SIGMA_E = 0.0233380;   % mass ratio sqrt(m_a/m_i) (correct = 0.0233380)
%% GRID PARAMETERS
NX      = 300;     % Realspace x-gridpoints
NY      = 300;     % Realspace y-gridpoints
LX      = 120;     % Size of the squared frequency domain
LY      = 120;     % Size of the squared frequency domain
NZ      = 1;      % number of perpendicular planes (parallel grid)
Q0      = 1.0;    % q factor ()
SHEAR   = 0.0;    % magnetic shear
EPS     = 0.0;    % inverse aspect ratio
P       = 4;
J       = 2;
%% TIME PARAMETERS
TMAX    = 1000;  % Maximal time unit
DT      = 1e-2;   % Time step
SPS0D   = 1/2;      % Sampling per time unit for profiler
SPS2D   = 1/2;      % Sampling per time unit for 2D arrays
SPS3D   = 1/2;      % Sampling per time unit for 3D arrays
SPS5D   = 1/50;  % Sampling per time unit for 5D arrays
JOB2LOAD= -1; % start from t=0 if <0, else restart from outputs_$job2load
%% OPTIONS AND NAMING
% Collision operator
% (LB:L.Bernstein, DG:Dougherty, SG:Sugama, LR: Lorentz, LD: Landau)
CO      = 'DG';
GKCO    = 1; % gyrokinetic operator
ABCO    = 1; % interspecies collisions
CLOS    = 0;   % Closure model (0: =0 truncation)
NL_CLOS = -1;   % nonlinear closure model (-2: nmax = jmax, -1: nmax = jmax-j, >=0 : nmax = NL_CLOS)
LINEARITY = 'nonlinear';   % activate non-linearity (is cancelled if KXEQ0 = 1)
% INIT options
INIT_ZF = 0; ZF_AMP = 0.0;
INIT_OPT     = 'phi'; 
ACT_ON_MODES = 'donothing';
%% OUTPUTS
W_DOUBLE = 0;
W_GAMMA  = 1; W_HF     = 1;
W_PHI    = 1; W_NA00   = 1;
W_DENS   = 1; W_TEMP   = 1;
W_NAPJ   = 1; W_SAPJ   = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% unused
KIN_E   = 1;         % Kinetic (1) or adiabatic (2) electron model
k_gB   = 1.0;       % Magnetic  gradient
k_cB   = 1.0;       % Magnetic  curvature
SG      = 0;         % Staggered z grids option
PMAXE   = P;    % Highest electron Hermite polynomial degree
JMAXE   = J;     % Highest ''       Laguerre ''
PMAXI   = P;     % Highest ion      Hermite polynomial degree
JMAXI   = J;     % Highest ''       Laguerre ''
KERN    = 0;   % Kernel model (0 : GK)
KX0KH   = 0; A0KH = 0; % Background phi mode to drive Ray-Tay inst.
KXEQ0   = 0;      % put kx = 0
KPAR    = 0.0;    % Parellel wave vector component
LAMBDAD = 0.0;
MU_X    = 1e-1; % Hyperdiffusivity coefficient
MU_Y    = 1e-1; % Hyperdiffusivity coefficient
NOISE0  = 1.0e-5;
BCKGD0  = 0.0;    % Init background
TAU     = 1.0;    % e/i temperature ratio
K_E     = 0.0;    % ES '''
MU_P    = 0.0;     % Hermite  hyperdiffusivity -mu_p*(d/dvpar)^4 f
MU_J    = 0.0;     % Laguerre hyperdiffusivity -mu_j*(d/dvperp)^4 f
% Compute processes distribution
Ntot = NP_P * NP_KX;
Nnodes = ceil(Ntot/48);
Nppn   = Ntot/Nnodes;
CLUSTER.NODES =  num2str(Nnodes);  % MPI process along p
CLUSTER.NTPN  =  num2str(Nppn); % MPI process along kx
CLUSTER.CPUPT = '1';        % CPU per task
%% Run file management scripts
setup
SBATCH_CMD = 'sbatch batch_script.sh\n';
write_sbash_marconi
if(mod(NP_P*NP_KX,48)~= 0)
    disp('WARNING : unused cores (ntot cores must be a 48 multiple)');
end
if(SUBMIT)
    [~,job_info_] = system('ssh ahoffman@login.marconi.cineca.it sh HeLaZ/wk/setup_and_run.sh');
    disp(job_info_);
    jobid_ = job_info_(21:27);
    if(CHAIN>0)
        for CHAIN_IDX = 1:CHAIN
            SBATCH_CMD = ['sbatch --dependency=afterok:',jobid_,' batch_script.sh\n'];
            disp(SBATCH_CMD);
            JOB2LOAD= JOB2LOAD+1;
%             DT = 1e-1;
%             NL_CLOS = -1;
            setup
            write_sbash_marconi
            [~,job_info_] = system('ssh ahoffman@login.marconi.cineca.it sh HeLaZ/wk/setup_and_run.sh');
            jobid_ = job_info_(21:27);
        end
    end
end
% system('rm fort*.90');
disp('done');
