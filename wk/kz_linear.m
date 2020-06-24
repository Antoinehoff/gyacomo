SIMID = 'kz_linear'; % Name of the simulations
addpath(genpath('../matlab')) % ... add 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% outputs options
OUTPUTS.nsave_0d = 0;
OUTPUTS.nsave_1d = 0;
OUTPUTS.nsave_2d = 100;
OUTPUTS.nsave_5d = 0;
OUTPUTS.write_Ni00    = '.true.';
OUTPUTS.write_moments = '.true.';
OUTPUTS.write_phi     = '.true.';
OUTPUTS.write_doubleprecision = '.true.';
OUTPUTS.resfile0      = '''results''';
%% Grid parameters
GRID.pmaxe = 81;
GRID.jmaxe = 20;
GRID.pmaxi = 81;
GRID.jmaxi = 20;
GRID.nkr   = 1;
GRID.krmin = 0.;
GRID.krmax = 0.;
GRID.nkz   = 10;
GRID.kzmin = 0.1/sqrt(2.0);
GRID.kzmax = 2.8/sqrt(2.0);
%% Model parameters
MODEL.nu      = 0.001;
MODEL.tau_e   = 1.0;
MODEL.tau_i   = 1.0;
MODEL.sigma_e = 0.0233380;
MODEL.sigma_i = 1.0;
MODEL.q_e     =-1.0;
MODEL.q_i     = 1.0;
MODEL.eta_n   = 1.0;
MODEL.eta_T   = 0.0;
MODEL.eta_B   = 0.5;
MODEL.lambdaD = 0.0;
%% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = '''RK4''';
BASIC.nrun              = 100000;
BASIC.dt                = 0.01;
BASIC.tmax              = 100;
INITIAL.initback_moments  = 0.01;
INITIAL.initnoise_moments = 0.;
INITIAL.iseed             = 42;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write input file
INPUT = write_fort90(OUTPUTS,GRID,MODEL,INITIAL,TIME_INTEGRATION,BASIC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run solver
nproc = 1;
EXEC  = ' ../bin/helaz ';
RUN   = ['mpirun -np ' num2str(nproc)];
CMD   = [RUN, EXEC, INPUT];
system(CMD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis and basic figures
SAVEFIG = 1;
helaz_analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%