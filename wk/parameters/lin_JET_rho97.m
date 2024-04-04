%% Parameters found in Parisi et al. 2020
% Jet shot 92174 parameters
BT0     = 1.9;      %[T]        Toroidal field @ 2.96m
Ip      = 1.4;      %[MA]       Plasma current @ 2.96m
PNBI    = 17.4;     %[MW]       NBI power
rhoi    = 0.27;     %[cm]       Ion gyroradius
R0      = 2.86;     %[m]        Major radius
a       = 0.91;     %[m]        F-T minor radius
Rc      = 2.91;     %[m]        ??
rc      = 0.89;     %[m]        ??
m_e     = 5.49e-4;  %[amu]      electron mass
m_i     = 2.014;    %[amu]      deuterium mass
% Dimless flux-tube parameters
nuee    = 0.83;     %[vti/a]    e-e collision frequ.
wTe     = 42;       %[a/L]      e-temp. gradient length
wTi     = 11;       %[a/L]      i-temp. gradient length
wNe     = 10;       %[a/L]      dens. gradient length
tau     = 1/0.56;   %Ti/Te      i-e temperature ratio
gE      = 0.56;     %[vti/a]    ExB shearing rate
roa     = 0.9743;   % r/a       Flux surface position
beta    = 0.0031;   % [8pi*ptot/B^2] with B = 1.99T
% Normalization       
% v0   = vth_i = sqrt(2*Ti/mi)                    
% rho0 = rho_i = vti/omegai
% Conversion factors from GYAC to paper results
freq_conv = a/R0 * sqrt(tau/2); % from R/c_s to a/vti
wave_conv = sqrt(2/tau);        % from rho_i  to rho_s
grad_conv = a/R0;               % from R/LT   to a/LT

%% Set simulation parameters
SIMID   = 'lin_JET_rho97';  % Name of the simulation
%% Set up physical parameters
NU      = 0.1;   % Not the true value 
TAU     = tau;         % i/e temperature ratio
K_Ne    = wNe/grad_conv;  % ele Density '''
K_Te    = wTe/grad_conv;  % ele Temperature '''
K_Ni    = wNe/grad_conv;  % ion Density gradient drive
K_Ti    = wTi/grad_conv;  % ion Temperature '''
SIGMA_E = sqrt(m_e/m_i);  % mass ratio sqrt(m_e/m_i) (e-H = 0.0233380)
NA      = 2;          % number of kinetic species
BETA    = beta;           % electron plasma beta
MHD_PD  = 1;
CO      = 'DG';       % Collision operator (LB:L.Bernstein, DG:Dougherty, SG:Sugama, LR: Lorentz, LD: Landau)
GKCO    = 1;          % Gyrokinetic operator
ABCO    = 1;          % INTERSPECIES collisions
COLL_KCUT= 100; % Cutoff for collision operator

%% GEOMETRY
% GEOMETRY= 's-alpha';
GEOMETRY= 'miller';
EPS     = a/R0;    % inverse aspect ratio
Q0      = 5.100;    % safety factor
SHEAR   = 3.360;    % magnetic shear
KAPPA   = 1.550;    % elongation
S_KAPPA = 0.949;
DELTA   = 0.263;    % triangularity
S_DELTA = 0.737;

%% Set up grid parameters
P = 4;
J = P/2;%P/2;
PMAX = P;                   % Hermite basis size
JMAX = J;                   % Laguerre basis size
NX = 16;                    % real space x-gridpoints
NY = 2;                     % real space y-gridpoints
LX = 2*pi/0.1;              % Size of the squared frequency domain in x direction
LY = 2*pi/0.2;             % Size of the squared frequency domain in y direction
NZ = 32;                    % number of perpendicular planes (parallel grid)
SG = 0;                     % Staggered z grids option
NEXC = 1;                   % To extend Lx if needed (Lx = Nexc/(kymin*shear))

%% TIME PARAMETERS
TMAX     = 15;  % Maximal time unit
DT       = 1e-3;   % Time step
DTSAVE0D = 0.5;      % Sampling time for 0D arrays
DTSAVE2D = -1;     % Sampling time for 2D arrays
DTSAVE3D = 0.5;      % Sampling time for 3D arrays
DTSAVE5D = 100;     % Sampling time for 5D arrays
JOB2LOAD = -1;     % Start a new simulation serie


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UNUSED PARAMETERS
% These parameters are usually not to play with in linear runs

%% OPTIONS
CLUSTER.TIME = '99:00:00';  % Allocation time hh:mm:ss
LINEARITY = 'linear';   % activate non-linearity (is cancelled if KXEQ0 = 1)
INIT_ZF   = 0;          % Initialize zero-field quantities
HRCY_CLOS = 'truncation';   % Closure model for higher order moments
DMAX      = -1;
NLIN_CLOS = 'truncation';   % Nonlinear closure model for higher order moments
NMAX      = 0;
KERN      = 0;   % Kernel model (0 : GK)
INIT_OPT  = 'phi';   % Start simulation with a noisy mom00/phi/allmom
NOISE0    = 1.0e-5; % Initial noise amplitude
BCKGD0    = 0.0e-5;    % Initial background
NUMERICAL_SCHEME = 'RK4'; % Numerical integration scheme (RK2,SSPx_RK2,RK3,SSP_RK3,SSPx_RK3,IMEX_SSP2,ARK2,RK4,DOPRI5)

%% OUTPUTS
W_DOUBLE = 1;     % Output flag for double moments
W_GAMMA  = 1;     % Output flag for gamma (Gyrokinetic Energy)
W_HF     = 1;     % Output flag for high-frequency potential energy
W_PHI    = 1;     % Output flag for potential
W_NA00   = 1;     % Output flag for nalpha00 (density of species alpha)
W_DENS   = 1;     % Output flag for total density
W_TEMP   = 1;     % Output flag for temperature
W_NAPJ   = 1;     % Output flag for nalphaparallel (parallel momentum of species alpha)
W_SAPJ   = 0;     % Output flag for saparallel (parallel current of species alpha)

%% Unused geometry
ZETA    = 0;    % squareness
S_ZETA  = 0;
PARALLEL_BC = 'dirichlet'; % Boundary condition for parallel direction ('dirichlet','periodic','shearless','disconnected')
SHIFT_Y  = 0.0;    % Shift in the periodic BC in z
NPOL     = 1;       % Number of poloidal turns
PB_PHASE = 0;

%% Diffusions
MU      = 0.0;    % Hyperdiffusivity coefficient
MU_X    = MU;     % Hyperdiffusivity coefficient in x direction
MU_Y    = MU;     % Hyperdiffusivity coefficient in y direction
N_HD    = 4;      % Degree of spatial-hyperdiffusivity
MU_Z    = 5.0;    % Hyperdiffusivity coefficient in z direction
HYP_V   = 'hypcoll'; % Kinetic-hyperdiffusivity model
MU_P    = 0.0;    % Hyperdiffusivity coefficient for Hermite
MU_J    = 0.0;    % Hyperdiffusivity coefficient for Laguerre
LAMBDAD = 0.0;    % Lambda Debye
K_gB   = 1.0;     % Magnetic gradient tuner  (1 is the real value)
K_cB   = 1.0;     % Magnetic curvature tuner (1 is the real value)
K_mB   = 1.0;     % mirror force tuner       (1 is the real value)
K_tB   = 1.0;     % trapping term tuner      (1 is the real value)
K_ldB  = 1.0;     % Landau damping tuner     (1 is the real value)
ADIAB_I = 0;          % adiabatic ion model
ADIAB_E = (NA==1);          % adiabatic electron model