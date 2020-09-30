%% Run MOLI for a time evolution of the moments at a given kperp
%% Move to MOLI workspace
cd ../../MoliSolver/MOLI/workspace/
%% Add paths
my_paths

%% Directories
ROOT = '/home/ahoffman/Documents/MoliSolver';
options.dirs.COSOlverdir = fullfile(ROOT,'COSOlver');
options.dirs.MOLIdir = fullfile(ROOT,'MOLI');

%% MOLI Physical and Main Parameters

% Solve DK AND/OR GK Linear Moment Hierarchy.
options.DKI  = 0;   % First-order DK
options.DKII = 0;   % Second-order DK
options.GK   = 1;   % Gyrokinetic GK
options.EM   = 0;   % Include Electromagnetic effects (only for GK=1)
options.GD   = 0;   % Gyro-Drift
options.GDI  = 0;   % First/second order

% Solve MOLI
options.MOLI = 1;    % 1 -> Solve MOLI, 0 -> off

% MOLI Solver
options.solver.solver = 3;

% MOLI Linear Fit Solver
options.LinFitSolver = 0;

%% Main parameter scan

% Closure by truncation
params.Pmaxi = GRID.pmaxi;           % parallel ion Hermite moments
params.Jmaxi = GRID.jmaxi;           % perpendicular ion Laguerre moments
params.Pmaxe = GRID.pmaxe;           % parallel electron Hermite moments
params.Jmaxe = GRID.jmaxe;           % perpendicular electron Laguerre moments

% w/wo ions
options.ions = 1;           % if ions are present -> 1, 0 otherwise

% Adiabatic electrons
options.electrons = 1;      % 0 ->  adiabatic electrons, 1 no adiabatic electrons

% w/wo soundwaves
options.sw = 1;             % 1 -> sound waves on, 0 -> put ion parallel velocity row/column to 0

%% Collision Operator Models and COSOlver Input Parameters
options.collI  = MODEL.CO;         % collI=-2 -> Dougherty, -1 -> COSOlver, 0 -> Lenard-Bernstein, other -> hyperviscosity exponent
options.collGK = 0;         % collDKGK =1 -> GK collision operator, else DK collision operator
options.COSOlver.GKE = 0;
options.COSOlver.GKI = 0;

% COSOlver Input Parameters (if collI = -1 only)
options.COSOlver.eecolls = 1;	   % 1 -> electron-electron collisions, 0 -> off
options.COSOlver.iicolls = 1;      % 1 -> ion-ion collisions, 0 -> off
options.COSOlver.eicolls = 1;      % 1 -> electron-ion collisions (e-i) on, 0 -> off
options.COSOlver.iecolls = 1;      % 1 -> ion-electron collisions (i-e) on, 0 -> off

% Collisional Coulomb sum bounds (only if collI = -1, i.e. Coulomb)
options.COSOlver.lmaxx = 10;                        % upper bound collision operator first sum first species
options.COSOlver.kmaxx = 10;                        % upper bound collision operator second sum first species
options.COSOlver.nmaxx = options.COSOlver.lmaxx;    % upper bound collision operator first sum second species
options.COSOlver.qmaxx = options.COSOlver.kmaxx;    % upper bound collision operator second sum second species

% Collsion FLR sum bounds
options.COSOlver.nemaxxFLR = 0;         % upper bound FLR electron collision
options.COSOlver.nimaxxFLR = 0;         % upper bound FLR ion collision

% Collision Operator Model
% Set electron/ion test and back-reaction model operator
%
%  0 => Coulomb Collisions
options.COSOlver.ETEST = 1; % 0 --> Buffer Operator, 1 --> Coulomb, 2 --> Lorentz
options.COSOlver.EBACK = 1;
options.COSOlver.ITEST = 1;
options.COSOlver.IBACK = 1;
options.COSOlver.ESELF = 1;
options.COSOlver.ISELF = 1;

options.COSOlver.OVERWRITE = 0;    % overwrite collisional matrices even if exists

options.COSOlver.cmd   = 'mpirun -np 6 ./bin/CO 2 2 2';
%% Physical Parameters

% Toroidal effects
options.magnetic        = 1;           % 1-> Add toroidal magnetic gradient drift resonance effects

% Physical Parameters
params.tau              = MODEL.tau_i; % Ti/Te
params.nu               = MODEL.nu;    % electron/ion collision frequency ... only for nu/ omega_pe < nuoveromegapemax (electron plasma frequency) [See Banks et al. (2017)]
params.nuoveromegapemax = inf;         % Maximum ratio between electron/ion collision frequency and electron plasma frequency [See Banks et al. (2017)]. Set to inf if not desired !!!
params.mu               = MODEL.sigma_e;   % sqrt(m_e/m_i)
params.kpar             = 0.0;         % normalized parallel wave number to the major radius
params.kperp            = kz_MOLI;  % normalized perpendicular toroidal wave number to the soundLarmor radius. Note: If ions ==0 (e.g. EPW), kperp --> b
params.kr               = kr_MOLI;  % Radial component of perpendicular vector
params.alphaD           = 0.0;         % (k*Debye length)^2
params.Rn               = MODEL.eta_n; % Major Radius / Background density gradient length
params.RTe              = MODEL.eta_T; % Major Radius * normalized kperp / Background electron temperature gradient length
params.RTi              = MODEL.eta_T; % Major Radius * normalized kperp / Background ion temperature gradient length
params.Rphi             = 0.0;         % Major Radius * normalized kperp / Background potentiel gradient length [presence of shear] - only for GK
params.betae            = 1e-6;        % Electron Beta plasma.

params.rhostar          = 1e-5;        % sound Larmor Radius/Major Radius ~ sqrt(Te)/(R_0*B).
params.n0               = INITIAL.initback_moments;        % initial density perturbation

params.gradB            = MODEL.eta_B;         % Magnetic field gradient
params.curvB            = MODEL.eta_B;         % Curvature of B
params.trappB           = 0.0;         % Trapping term

%% MOLI Options

% Save data in results dir
options.save    = 0;
options.verbose = 0;
options.dbg     = 0;

options.DR        = 0;      % 1 -> Solve kinetic dispersion relation,
options.KineticDR = 0;      % Solve kinetic dispersion relation (Landau integral) for the given theory

% Compute the kinetic susceptibility for EPW only
options.SPTBLTY = 0;

options.nharm   = 1;        % Number of harmonics in disp. rel. 1 and 4

wlim = 5.0;
options.DRsolver.wr_min = -wlim;     % Minimum of real part.
options.DRsolver.wr_max =  wlim;     % Maximum of real part.
options.DRsolver.wi_min = -wlim;     % Minimum of imag part.
options.DRsolver.wi_max =  wlim;     % Maximum of imag part.
options.DRsolver.nw     =  300;      % Grid resolution

% Disp. Rel. Options
options.FLRmodel    = 0;   % 1 -> Truncated Laguerre, 0 -> Exact representation
options.FluidLandau = 0;   % 1 -> Add Landau Fluid Closure to Fluid Dispersion Relation, 0 -> off
options.deltaLandau = 0;   % 1 -> Hammet-Perkins closure on, 0 -> off

% Fluid dispersion relation
options.FluidDR = 0;	     % Solve annamaria's fluid equations
options.Fluid.sITGmm = 0;

% Define scan parameters
options.fscan = 0;         % 1 -> peform scan over scan.list, 0-> off

options.scan.list = {};% List of scan parameters. If empty, solve MOLI with params

% Time-Evolution Problem [Solver==3] ...
options.solver.TimeSolver.dt         = BASIC.dt;		% timestep of time evolution (R/c_s or 1/(k v/the) units)
options.solver.TimeSolver.tmax       = BASIC.tmax;
options.solver.TimeSolver.Trun       = BASIC.tmax;		  % total time to run time evolution
options.solver.TimeSolver.t_fit_min  = 0.05;    % Phase-Mixing fit Lower time limit
options.solver.TimeSolver.t_fit_max  = 8;       % Phase-Mixing fit Upper time limit
options.solver.TimeSolver.en_fit_min = 0.15;    % Entropy Mode fit Lower time limit
options.solver.TimeSolver.en_fit_max = 0.3;     % Entropy Mode fit Upper time limit
options.solver.TimeSolver.movie      = 0;       % Display movie if 1, last frame otherwise
options.solver.TimeSolver.save       = 0;       % 1 --> save during fscan, Warning: memory storage


%% Run MOLI
% Solve the MOLI
[results,params,options] = MOLI_Control(params,options);

%% Return to HeLaZ workspace
cd ../../../HeLaZ/wk
