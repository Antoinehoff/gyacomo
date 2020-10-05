clear all; close all;
addpath(genpath('../matlab')) % ... add
default_plots_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 1e-2;   % Collision frequency
TAU     = 1.0;    % e/i temperature ratio
ETAB    = 0.1;    % Magnetic gradient
ETAN    = 1.0;    % Density gradient
ETAT    = 0.0;    % Temperature gradient
MU      = 1e-2;   % Hyper diffusivity coefficient
LAMBDAD = 0.0;
%% GRID PARAMETERS
N       = 64;     % Frequency gridpoints (Nr = N/2)
L       = 10;      % Size of the squared frequency domain
PMAXE   = 05;     % Highest electron Hermite polynomial degree
JMAXE   = 02;     % Highest ''       Laguerre ''
PMAXI   = 05;     % Highest ion      Hermite polynomial degree
JMAXI   = 02;     % Highest ''       Laguerre ''
KPAR    = 0.0;    % Parellel wave vector component
%% TIME PARAMETERS 
TMAX    = 100.0;    % Maximal time unit
DT      = 1e-4;   % Time step
SPS     = 10;     % Sampling per time unit
RESTART = 1;      % To restart from last checkpoint
JOB2LOAD= 0;
%% OPTIONS
SIMID   = 'gvskr';  % Name of the simulation
NON_LIN = 0;   % activate non-linearity (is cancelled if KREQ0 = 1)
CO      = -2;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
DK      = 0;   % Drift kinetic model (put every kernel to 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ________________________________________________________________________
% Naming and creating input file
params   = ['P_',num2str(PMAXE),'_J_',num2str(JMAXE),...
    '_nB_',num2str(ETAB),'_nN_',num2str(ETAN),'_mu_%0.0e_'];
params   = sprintf(params,MU);
if ~NON_LIN; params = ['lin_',params]; end;
if  DK;      params = ['DK_',params]; end;
resolution = ['_',num2str(N/2),'x',num2str(N),'_'];
gridname   = ['L_',num2str(L),'_'];
BASIC.SIMID = [SIMID,resolution,gridname,params];
BASIC.nrun       = 1e8;
BASIC.dt         = DT;   
BASIC.tmax       = TMAX;    %time normalized to 1/omega_pe
if RESTART; BASIC.RESTART = '.true.'; else; BASIC.RESTART = '.false.';end;
OUTPUTS.nsave_0d = 0;
OUTPUTS.nsave_1d = 0;
OUTPUTS.nsave_2d = floor(1.0/SPS/DT);
OUTPUTS.nsave_5d = floor(1.0/SPS/DT);
OUTPUTS.write_Ni00    = '.true.';
OUTPUTS.write_moments = '.true.';
OUTPUTS.write_phi     = '.true.';
OUTPUTS.write_non_lin = '.true.';
if NON_LIN == 0; OUTPUTS.write_non_lin = '.false.'; end;
OUTPUTS.write_doubleprecision = '.true.';
OUTPUTS.resfile0      = ['''',BASIC.SIMID,''''];
OUTPUTS.rstfile0      = ['''','../checkpoint/cp_',BASIC.SIMID,''''];
OUTPUTS.job2load      = JOB2LOAD;
% Grid parameters
GRID.pmaxe = PMAXE;  % Electron Hermite moments
GRID.jmaxe = JMAXE;  % Electron Laguerre moments 
GRID.pmaxi = PMAXI;  % Ion Hermite moments
GRID.jmaxi = JMAXI;  % Ion Laguerre moments
GRID.Nr    = N; % r grid resolution
GRID.Lr    = L; % r length
GRID.Nz    = N; % z ''
GRID.Lz    = L; % z ''
GRID.kpar  = KPAR;
% Model parameters
MODEL.CO      = CO;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
if DK;      MODEL.DK      = '.true.'; else; MODEL.DK      = '.false.';end;
if NON_LIN; MODEL.NON_LIN = '.true.'; else; MODEL.NON_LIN = '.false.';end;
MODEL.mu      = MU;
MODEL.nu      = NU; % hyper diffusive coefficient nu for HW
% temperature ratio T_a/T_e
MODEL.tau_e   = TAU;
MODEL.tau_i   = TAU;
% mass ratio sqrt(m_a/m_i)
MODEL.sigma_e = 0.0233380;
MODEL.sigma_i = 1.0;
% charge q_a/e
MODEL.q_e     =-1.0;
MODEL.q_i     = 1.0;
% gradients L_perp/L_x
MODEL.eta_n   = ETAN;        % source term kappa for HW
MODEL.eta_T   = ETAT;        % Temperature
MODEL.eta_B   = ETAB;        % Magnetic
MODEL.lambdaD = LAMBDAD;
% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = '''RK4''';
INITIAL.only_Na00         = '.false.';
INITIAL.initback_moments  = 1.0e-4;
INITIAL.initnoise_moments = 5.0e-5;
INITIAL.iseed             = 42;
INITIAL.selfmat_file = ...
    ['''../iCa/self_Coll_GKE_0_GKI_0_ESELF_1_ISELF_1_Pmaxe_',num2str(GRID.pmaxe),...
    '_Jmaxe_',num2str(GRID.jmaxe),'_Pmaxi_',num2str(GRID.pmaxi),'_Jmaxi_',...
    num2str(GRID.jmaxi),'_pamaxx_10.h5'''];
INITIAL.eimat_file = ...
    ['''../iCa/ei_Coll_GKE_0_GKI_0_ETEST_1_EBACK_1_Pmaxe_',num2str(GRID.pmaxe),...
    '_Jmaxe_',num2str(GRID.jmaxe),'_Pmaxi_',num2str(GRID.pmaxi),'_Jmaxi_',...
    num2str(GRID.jmaxi),'_pamaxx_10_tau_1.0000_mu_0.0233.h5'''];
INITIAL.iemat_file = ...
    ['''../iCa/ie_Coll_GKE_0_GKI_0_ITEST_1_IBACK_1_Pmaxe_',num2str(GRID.pmaxe),...
    '_Jmaxe_',num2str(GRID.jmaxe),'_Pmaxi_',num2str(GRID.pmaxi),'_Jmaxi_',...
    num2str(GRID.jmaxi),'_pamaxx_10_tau_1.0000_mu_0.0233.h5'''];

%% Compile and write input file
INPUT = write_fort90(OUTPUTS,GRID,MODEL,INITIAL,TIME_INTEGRATION,BASIC);
nproc = 1;
MAKE  = 'cd ..; make; cd wk';
system(MAKE);
%%
disp(['Set up ', BASIC.SIMID]);