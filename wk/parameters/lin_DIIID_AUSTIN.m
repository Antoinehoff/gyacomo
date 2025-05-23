%% Reference values
prof_folder = ['parameters/profiles/DIIID_Austin_et_al_2019/',TRIANG,'/'];
% Get geometry parameters
GEOMETRY= 'miller';
if READPROF
    geom    = get_miller_GENE_py(prof_folder,rho);
else
    geom.q0     = 3.25; geom.shat   = 3.65;
    geom.kappa  = 1.66; geom.s_kappa= 0.77;
    geom.delta  = 0.34; geom.s_delta= 1.13;
    if strcmp(TRIANG,'NT')
        geom.delta = -geom.delta;
    end
    geom.zeta   =-0.038;geom.s_zeta =-0.13;
    geom.trpeps = 1/3.15;
    geom.Bref   = 1;
    geom.R0     = 3.15;
    geom.a      = 1;
    geom.Lref   = geom.R0;
end
Q0      = geom.q0;    % safety factor
SHEAR   = geom.shat;    % magnetic shear
EPS     = geom.trpeps;    % inverse aspect ratio
KAPPA   = geom.kappa;     % elongation
S_KAPPA = geom.s_kappa;
DELTA   = geom.delta;    % triangularity
S_DELTA = geom.s_delta;
ZETA    = geom.zeta;    % squareness
S_ZETA  = geom.s_zeta;
PARALLEL_BC = 'dirichlet'; % Boundary condition for parallel direction ('dirichlet','periodic','shearless','disconnected')
% See Austin et al. 2019, negative triangularity exp. in DIIID
q_e= 1.602e-19;
Z  = 2;
m_i= Z*1.67e-27;
m_e= 9.11e-31;
BT = geom.Bref; %T
Ip = 0.9; %MA
R0 = geom.R0; %m
a0 = geom.a;%m
P_in = 13;  %input power in MW
Lref = geom.Lref;
Bref = abs(geom.Bref);
mref = m_i;
% Get values at a given location :
[params,profiles] = get_param_from_profiles(prof_folder,rho,Lref,mref,Bref,READPROF);
K_Ne   = params.K_Ne;
K_Ni   = params.K_Ne;
K_Te   = params.K_Te;
K_Ti   = params.K_Ti;
TAU    = params.TAU;
NU     = params.NU;
BETA   = params.BETA;
SIGMA_E= sqrt(m_e/mref);
EXBRATE= params.EXBRATE;
SHIFT_Y = 0.0;    % Shift in the periodic BC in z
NPOL   = 1;       % Number of poloidal turns
PB_PHASE = 0;
% Name of the simulation folder
SIMID   = ['lin_DIIID_Austin',TRIANG,'_',num2str(rho)]; 
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
DTSAVE5D = TMAX/2;     % Sampling time for 5D arrays
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
NOISE0  = 1.0e-5; % Initial noise amplitude
BCKGD0  = 0.0e-5;    % Initial background
K_gB   = 1.0;     % Magnetic gradient tuner  (1 is the real value)
K_cB   = 1.0;     % Magnetic curvature tuner (1 is the real value)
K_mB   = 1.0;     % mirror force tuner       (1 is the real value)
K_tB   = 1.0;     % trapping term tuner      (1 is the real value)
K_ldB  = 1.0;     % Landau damping tuner     (1 is the real value)
COLL_KCUT = 1; % Cutoff for collision operator
ADIAB_I = 0;          % adiabatic ion model