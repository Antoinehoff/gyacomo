% Taken from Table 2, Kennedy et al. 2023 preprint on STEP linear
% simulations. rho = 0.7
%% Reference values
% Get geometry parameters
GEOMETRY= 'miller';
Q0      = 4.0;    % safety factor
SHEAR   = 1.56;    % magnetic shear
EPS     = 0.56;    % inverse aspect ratio
KAPPA   = 2.57;     % elongation
S_KAPPA = 0.19;
DELTA   = 0.32;    % triangularity
S_DELTA = 0.54;
ZETA    = 0;    % squareness
S_ZETA  = 0;
PARALLEL_BC = 'dirichlet'; % Boundary condition for parallel direction ('dirichlet','periodic','shearless','disconnected')
% See Austin et al. 2019, negative triangularity exp. in DIIID
q_e= 1.602e-19;
Z  = 2;
m_i= Z*1.67e-27;
m_e= 9.11e-31;
BT = 3.2; %T
Ip = 20.9; %MA
R0 = 3.6; %m
a0 = 2;%m
P_in = 300;  %input power in MW
Lref = R0;
Bref = abs(BT);
mref = m_i;
T_e  = 18.0; % KeV
n_e  = 20.05; %10^19 m-3
% Get values at a given location :
K_Ne   = 1.54;
K_Ni   = K_Ne;
K_Te   = 1.77;
K_Ti   = 1.96;
TAU    = 1.0;
nuGENE = 2.3031E-5*Lref*(n_e)/(T_e)^2*(24.-log(sqrt(n_e*1.0E13)/T_e*0.001));
NU     = 0.05;%3/8*sqrt(pi/2)*sqrt(TAU)*nuGENE;
BETA   = 0.07;%403.0E-5*n_e*T_e/(Bref*Bref);
SIGMA_E= sqrt(m_e/mref);
EXBRATE= 0;
SHIFT_Y = 0.0;    % Shift in the periodic BC in z
NPOL   = 1;       % Number of poloidal turns
PB_PHASE = 0;
% Name of the simulation folder
SIMID   = 'lin_STEP_EC_HD_psi71'; 
%% Set up physical parameters
CLUSTER.TIME = '99:00:00';  % Allocation time hh:mm:ss
NA      = 2;          % number of kinetic species
ADIAB_E = (NA==1);          % adiabatic electron model
MHD_PD  = 1;
%% Set up grid parameters
P = 2;
J = P/2;%P/2;
PMAX = P;                   % Hermite basis size
JMAX = J;                   % Laguerre basis size
NX = 6;                    % real space x-gridpoints
NY = 2;                     % real space y-gridpoints
LX = 2*pi/0.1;              % Size of the squared frequency domain in x direction
LY = 2*pi/0.3;             % Size of the squared frequency domain in y direction
NZ = 32;                    % number of perpendicular planes (parallel grid)
SG = 0;                     % Staggered z grids option
NEXC = 1;                   % To extend Lx if needed (Lx = Nexc/(kymin*shear))
%% TIME PARAMETERS
TMAX     = 15;  % Maximal time unit
DT       = 1e-4;   % Time step
DTSAVE0D = 0.5;      % Sampling time for 0D arrays
DTSAVE2D = -1;     % Sampling time for 2D arrays
DTSAVE3D = 0.5;      % Sampling time for 3D arrays
DTSAVE5D = 100;     % Sampling time for 5D arrays
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
INIT_OPT  = 'phi';   % Start simulation with a noisy mom00/phi/allmom
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
NOISE0  = 1.0e-5; % Initial noise amplitude
BCKGD0  = 0.0e-5;    % Initial background
K_gB   = 1.0;     % Magnetic gradient tuner  (1 is the real value)
K_cB   = 1.0;     % Magnetic curvature tuner (1 is the real value)
K_mB   = 1.0;     % mirror force tuner       (1 is the real value)
K_tB   = 1.0;     % trapping term tuner      (1 is the real value)
K_ldB  = 1.0;     % Landau damping tuner     (1 is the real value)
COLL_KCUT = 1; % Cutoff for collision operator
ADIAB_I = 0;          % adiabatic ion model