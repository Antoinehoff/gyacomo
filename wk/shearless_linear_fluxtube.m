RUN = 1; % To run or just to load
addpath(genpath('../matlab')) % ... add
default_plots_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
CLUSTER.TIME  = '99:00:00'; % allocation time hh:mm:ss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 0.0;       % Collision frequency
TAU     = 1.0;       % e/i temperature ratio
K_N     = 2.22;      % Density gradient drive
K_T     = 6.0;       % Temperature '''
SIGMA_E = 0.05196;   % mass ratio sqrt(m_a/m_i) (correct = 0.0233380)
KIN_E   = 1;         % Kinetic (1) or adiabatic (0) electron model
%% GRID PARAMETERS
NX      = 1;         % real space x-gridpoints
NY      = 2;         %     ''     y-gridpoints
LX      = 0;         % Size of the squared frequency domain
LY      = 2*pi/0.25; % Size of the squared frequency domain
NZ      = 24;        % number of perpendicular planes (parallel grid)
Q0      = 2.7;       % safety factor
SHEAR   = 0.0;       % magnetic shear
EPS     = 0.18;      % inverse aspect ratio
%% TIME PARMETERS
TMAX    = 10;  % Maximal time unit
DT      = 1e-3;   % Time step
SPS0D   = 1;      % Sampling per time unit for 2D arrays
SPS2D   = 0;      % Sampling per time unit for 2D arrays
SPS3D   = 10;      % Sampling per time unit for 2D arrays
SPS5D   = 1/100;    % Sampling per time unit for 5D arrays
SPSCP   = 0;    % Sampling per time unit for checkpoints
JOB2LOAD= -1;
%% OPTIONS
SIMID   = 'shearless_fluxtube';  % Name of the simulation
% Collision operator
% (0 : L.Bernstein, 1 : Dougherty, 2: Sugama, 3 : Pitch angle, 4 : Full Couloumb ; +/- for GK/DK)
CO      = 1;
INIT_ZF = 0; ZF_AMP = 0.0;
CLOS    = 0;   % Closure model (0: =0 truncation, 1: gyrofluid closure (p+2j<=Pmax))
NL_CLOS =-1;   % nonlinear closure model (-2:nmax=jmax; -1:nmax=jmax-j; >=0:nmax=NL_CLOS)
KERN    = 0;   % Kernel model (0 : GK)
%% OUTPUTS
W_DOUBLE = 0;
W_GAMMA  = 1; W_HF     = 1;
W_PHI    = 1; W_NA00   = 1;
W_DENS   = 1; W_TEMP   = 1;
W_NAPJ   = 1; W_SAPJ   = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unused
HD_CO   = 0.5;    % Hyper diffusivity cutoff ratio
kmax    = Nx*pi/Lx;% Highest fourier mode
NU_HYP  = 0.0;    % Hyperdiffusivity coefficient
MU      = NU_HYP/(HD_CO*kmax)^4; % Hyperdiffusivity coefficient
MU_P    = 0.0;     % Hermite  hyperdiffusivity -mu_p*(d/dvpar)^4 f
MU_J    = 0.0;     % Laguerre hyperdiffusivity -mu_j*(d/dvperp)^4 f
K_E     = 0.00;   % Electrostat '''
GRADB   = 1.0;
CURVB   = 1.0;
INIT_BLOB = 0; WIPE_TURB = 0; WIPE_ZF = 0; INIT_PHI= 0;
NOISE0  = 0.0e-4; % Init noise amplitude
BCKGD0  = 1.0e-4; % Init background
LAMBDAD = 0.0;
KXEQ0   = 0;      % put kx = 0
NON_LIN = 0;   % activate non-linearity (is cancelled if KXEQ0 = 1)
%% PARAMETER SCANS

if 1
%% Parameter scan over PJ
% PA = [2 4];
% JA = [1 2];
PA = [4];
JA = [2];
DTA= DT*ones(size(JA));%./sqrt(JA);
% DTA= DT;
mup_ = MU_P;
muj_ = MU_J;
Nparam = numel(PA);
param_name = 'PJ';
gamma_Ni00 = zeros(Nparam,floor(Nx/2)+1);
gamma_Nipj = zeros(Nparam,floor(Nx/2)+1);
gamma_phi  = zeros(Nparam,floor(Nx/2)+1);
% Ni00_ST  = zeros(Nparam,floor(Nx/2)+1,TMAX/SPS3D);
%  PHI_ST  = zeros(Nparam,floor(Nx/2)+1,TMAX/SPS3D);
for i = 1:Nparam
    % Change scan parameter
    PMAXE = PA(i); PMAXI = PA(i);
    JMAXE = JA(i); JMAXI = JA(i);
    DT = DTA(i);
    setup
%     system(['rm fort*.90']);
    % Run linear simulation
    if RUN
        system(['cd ../results/',SIMID,'/',PARAMS,'/; ./../../../bin/helaz 0; cd ../../../wk'])
    end
%     Load and process results
    %%
    filename = ['../results/',SIMID,'/',PARAMS,'/outputs_00.h5'];
    load_results
end

end

if 0
%% Trajectories of some modes
figure;
for i = 1:10:Nx/2+1
    semilogy(Ts3D,squeeze(abs(Ne00(i,2,1,:))),'DisplayName',['k=',num2str(kx(i))]); hold on;
end
end