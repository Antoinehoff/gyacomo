%clear all;
addpath(genpath('../matlab')) % ... add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLUSTER PARAMETERS
CLUSTER.TIME  = '01:00:00';   % allocation time hh:mm:ss
CLUSTER.NODES = '1';        % MPI process
CLUSTER.CPUPT = '1';        % CPU per task
CLUSTER.NTPN  = '4';       % N tasks per node
CLUSTER.MEM   = '2000';       % Memory
CLUSTER.JNAME = 'test_HeLaZ'; % Job name
USERNAME      = 'ahoffman';   % username at ppb110 for folder naming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 1.0;   % Collision frequency
ETAB    = 0.5;   % Magnetic gradient
NU_HYP  = 0.1;   % Hyperdiffusivity coefficient
%% GRID PARAMETERS
N       = 150;     % Frequency gridpoints (Nkr = N/2)
L       = 70;     % Size of the squared frequency domain
P       = 2;       % Electron and Ion highest Hermite polynomial degree
J       = 1;       % Electron and Ion highest Laguerre polynomial degree
MU_P    = 0;     % Hermite  hyperdiffusivity -mu_p*(d/dvpar)^4 f
MU_J    = 0;     % Laguerre hyperdiffusivity -mu_j*(d/dvperp)^4 f
%% TIME PARAMETERS
TMAX    = 150;  % Maximal time unit
DT      = 2e-2;  % Time step
SPS0D   = 1;      % Sampling per time unit for profiler
SPS2D   = 1/10;   % Sampling per time unit for 2D arrays
SPS5D   = 1/10;  % Sampling per time unit for 5D arrays
SPSCP   = 0;     % Sampling per time unit for checkpoints
RESTART = 0;     % To restart from last checkpoint
JOB2LOAD= 0;
%% OPTIONS
SIMID   = 'test_low_sampling';  % Name of the simulation
SIMID   = sprintf(SIMID,NU);
CO      = -3;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty, -3 : GK Dougherty)
CLOS    = 0;   % Closure model (0: =0 truncation, 1: semi coll, 2: Copy closure J+1 = J, P+2 = P)
KERN    = 0;   % Kernel model (0 : GK)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fixed parameters (for current study)
KR0KH   = 0; A0KH = 0; % Background phi mode to drive Ray-Tay inst. (not implemented)
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
setup % Write the input script "fort.90" with desired parameters

%% Write the sh script to launch the job on slurm for PPB110
INPUT = 'setup_and_run.sh';
fid = fopen(INPUT,'wt');
SCRATCH_SIMDIR = ['/scratch/',USERNAME,'/HeLaZ']; % Path to the simulation directory in the scratch
% Writing of the script
fprintf(fid,[...
'#!/bin/bash\n',...
'mkdir -p ',SCRATCH_SIMDIR,'/wk\n',...
...
'cd ',SCRATCH_SIMDIR,'/wk\n',...
...
'mkdir -p ', BASIC.RESDIR,'\n',...
'cd ',BASIC.RESDIR,'\n',...
'cp $HOME/HeLaZ/wk/fort.90 .\n',...
'cp $HOME/HeLaZ/wk/batch_script.sh .\n',...
...
'sbatch batch_script.sh\n',...
'echo $',SCRATCH_SIMDIR,'/results/',BASIC.SIMID,'/',BASIC.PARAMS,'/out.txt']);

fclose(fid);
system(['cp setup_and_run.sh ',BASIC.RESDIR,'/.']);

%% Write the sbatch script for PPB110
INPUT = 'batch_script.sh';
fid = fopen(INPUT,'wt');

fprintf(fid,[...
'#!/bin/bash\n',...
'#SBATCH --job-name=',CLUSTER.JNAME,'\n',...
'#SBATCH --time=', CLUSTER.TIME,'\n',...
'#SBATCH --nodes=', CLUSTER.NODES,'\n',...
'#SBATCH --cpus-per-task=', CLUSTER.CPUPT,'\n',...
'#SBATCH --ntasks-per-node=', CLUSTER.NTPN,'\n',...
'#SBATCH --mem-per-cpu=', CLUSTER.MEM,'\n',...
'#SBATCH --error=err.txt\n',...
'#SBATCH --output=out.txt\n',...
'module purge\n',...
'module load PrgEnv-intel/17.0\n',...
'srun  ./../../../bin/helaz']);

fclose(fid);
system(['cp batch_script.sh ',BASIC.RESDIR,'/.']);

disp('done');
