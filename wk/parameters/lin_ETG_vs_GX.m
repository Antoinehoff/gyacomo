%% Set simulation parameters
SIMID = 'lin_ETG_vs_GX'; % Name of the simulation

%% Set up physical parameters
CLUSTER.TIME = '99:00:00';  % Allocation time hh:mm:ss
NU      = 0.005;            % Collision frequency
TAU     = 1.0;              % e/i temperature ratio
K_Ne    = 3.0;           % ele Density
K_Te    = 9.0;           % ele Temperature
K_Ni    = 3.0;             % ion Density gradient drive
K_Ti    = 9.0;             % ion Temperature
SIGMA_E = 0.0165;        % mass ratio sqrt(m_a/m_i) (correct e-H = 0.0233380)
NA      = 2;                % number of kinetic species
ADIAB_E = (NA==1);          % adiabatic electron model
BETA    = 0.05;             % electron plasma beta
EXBRATE = 0.0;              % Background ExB shear flow
%% Set up grid parameters
PMAX = 32;                   % Hermite basis size
JMAX = 24;                   % Laguerre basis size
NX = 12;                     % real space x-gridpoints
NY =  2;                    % real space y-gridpoints
LX = 1.8;                % Size of the squared frequency domain in x direction
LY = 2*pi/10.0;              % Size of the squared frequency domain in y direction
NZ = 64;                    % number of perpendicular planes (parallel grid)
SG = 0;                     % Staggered z grids option
NEXC = 0;                   % To extend Lx if needed (Lx = Nexc/(kymin*shear))

%% GEOMETRY
GEOMETRY= 's-alpha';
% GEOMETRY= 'miller';
EPS     = 0.16667;   % inverse aspect ratio
Q0      = 2.0;    % safety factor
SHEAR   = 1.0;    % magnetic shear
KAPPA   = 1.0;    % elongation
DELTA   = 0.0;    % triangularity
ZETA    = 0.0;    % squareness
PARALLEL_BC = 'dirichlet'; % Boundary condition for parallel direction ('dirichlet','periodic','shearless','disconnected')
SHIFT_Y = 0.0;    % Shift in the periodic BC in z
NPOL   = 1;       % Number of poloidal turns

%% TIME PARAMETERS
TMAX     = 1.0;  % Maximal time unit
DT       = 1e-4;   % Time step
DTSAVE0D = 0.01;      % Sampling time for 0D arrays
DTSAVE2D = -1;     % Sampling time for 2D arrays
DTSAVE3D = 0.01;      % Sampling time for 3D arrays
DTSAVE5D = 0.01;     % Sampling time for 5D arrays
JOB2LOAD = -1;     % Start a new simulation serie

%% OPTIONS
LINEARITY = 'linear';   % activate non-linearity (is cancelled if KXEQ0 = 1)
CO        = 'DG';       % Collision operator (LB:L.Bernstein, DG:Dougherty, SG:Sugama, LR: Lorentz, LD: Landau)
GKCO      = 1;          % Gyrokinetic operator
ABCO      = 1;          % INTERSPECIES collisions
INIT_ZF   = 0;          % Initialize zero-field quantities
HRCY_CLOS = 'truncation';   % Closure model for higher order moments
DMAX      = -1;
NLIN_CLOS = 'truncation';   % Nonlinear closure model for higher order moments
NMAX      = 0;
KERN      = 0;   % Kernel model (0 : GK)
INIT_OPT  = 'mom00';   % Start simulation with a noisy mom00/phi/allmom
NUMERICAL_SCHEME = 'RK4'; % Numerical integration scheme (RK2,SSPx_RK2,RK3,SSP_RK3,SSPx_RK3,IMEX_SSP2,ARK2,RK4,DOPRI5)

%% OUTPUTS
W_DOUBLE = 1;     % Output flag for double moments
W_GAMMA  = 1;     % Output flag for gamma (Gyrokinetic Energy)
W_HF     = 1;     % Output flag for high-frequency potential energy
W_PHI    = 1;     % Output flag for potential
W_NA00   = 1;     % Output flag for nalpha00 (density of species alpha)
W_DENS   = 1;     % Output flag for total density
W_FVEL   = 1;     % Output flag for fluid velocity
W_TEMP   = 1;     % Output flag for temperature
W_NAPJ   = 1;     % Output flag for nalphaparallel (parallel momentum of species alpha)
W_SAPJ   = 0;     % Output flag for saparallel (parallel current of species alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UNUSED PARAMETERS
% These parameters are usually not to play with in linear runs
MU      = 0.0;    % Hyperdiffusivity coefficient
MU_X    = MU;     % Hyperdiffusivity coefficient in x direction
MU_Y    = MU;     % Hyperdiffusivity coefficient in y direction
N_HD    = 4;      % Degree of spatial-hyperdiffusivity
MU_Z    = 5.0;    % Hyperdiffusivity coefficient in z direction
HYP_V   = 'hypcoll'; % Kinetic-hyperdiffusivity model
MU_P    = 0.0;    % Hyperdiffusivity coefficient for Hermite
MU_J    = 0.0;    % Hyperdiffusivity coefficient for Laguerre
LAMBDAD = 0.0;    % Lambda Debye
NOISE0  = 0.0e-5; % Initial noise amplitude
BCKGD0  = 1.0e-5;    % Initial background
K_gB   = 1.0;     % Magnetic gradient tuner  (1 is the real value)
K_cB   = 1.0;     % Magnetic curvature tuner (1 is the real value)
K_mB   = 1.0;     % mirror force tuner       (1 is the real value)
K_tB   = 1.0;     % trapping term tuner      (1 is the real value)
K_ldB  = 1.0;     % Landau damping tuner     (1 is the real value)
COLL_KCUT = 1; % Cutoff for collision operator
S_KAPPA = 0.0;
S_DELTA = 0.0;
S_ZETA  = 0.0;
PB_PHASE= 0;
ADIAB_I = 0;
MHD_PD  = 0;