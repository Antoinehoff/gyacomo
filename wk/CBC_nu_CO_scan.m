% NU_a = [0.05 0.15 0.25 0.35 0.45];
NU_a = [0:0.1:0.5];
g_max= NU_a*0;
g_avg= NU_a*0;
g_std= NU_a*0;
k_max= NU_a*0;
CO      = 'DG';

K_T   = 7;
DT    = 5e-3;
TMAX  = 20;
ky_   = 0.3;
SIMID = 'linear_CBC_nu_scan_kT_7_ky_0.3_DGGK';  % Name of the simulation
% SIMID = 'linear_CBC_nu_scan_kT_11_ky_0.3_DGGK';  % Name of the simulation
RUN   = 1;
figure

for P = [6 8 10]

i=1;
for NU = NU_a
    
%Set Up parameters
for j = 1
    CLUSTER.TIME  = '99:00:00'; % allocation time hh:mm:ss
    TAU     = 1.0;    % e/i temperature ratio
    K_N = 2.22; K_Ne = K_N;
    K_Te     = K_T;            % Temperature '''
    SIGMA_E = 0.0233380;   % mass ratio sqrt(m_a/m_i) (correct = 0.0233380)
    KIN_E   = 0;     % 1: kinetic electrons, 2: adiabatic electrons
    BETA    = 0e-1;     % electron plasma beta 
    J = P/2;
    PMAXE   = P; JMAXE   = J;
    PMAXI   = P; JMAXI   = J;
    NX      = 12;    % real space x-gridpoints
    NY      = 2;     %     ''     y-gridpoints
    LX      = 2*pi/0.1;   % Size of the squared frequency domain
    LY      = 2*pi/ky_;
    NZ      = 16;    % number of perpendicular planes (parallel grid)
    NPOL    = 1; SG = 0;
    GEOMETRY= 's-alpha';
    Q0      = 1.4;    % safety factor
    SHEAR   = 0.8;    % magnetic shear (Not implemented yet)
    EPS     = 0.18;    % inverse aspect ratio
    SPS0D = 1; SPS2D = 0; SPS3D = 1;SPS5D= 1/5; SPSCP = 0;
    JOB2LOAD= -1;
    LINEARITY = 'linear';   % activate non-linearity (is cancelled if KXEQ0 = 1)
    GKCO    = 1; % gyrokinetic operator
    ABCO    = 1; % interspecies collisions
    INIT_ZF = 0; ZF_AMP = 0.0;
    CLOS    = 0;   % Closure model (0: =0 truncation, 1: v^Nmax closure (p+2j<=Pmax))s
    NL_CLOS = 0;   % nonlinear closure model (-2:nmax=jmax; -1:nmax=jmax-j; >=0:nmax=NL_CLOS)
    KERN    = 0;   % Kernel model (0 : GK)
    INIT_OPT= 'phi';   % Start simulation with a noisy mom00/phi/allmom
    W_DOUBLE = 1;
    W_GAMMA  = 1; W_HF     = 1;
    W_PHI    = 1; W_NA00   = 1;
    W_DENS   = 1; W_TEMP   = 1;
    W_NAPJ   = 1; W_SAPJ   = 0;
    HD_CO   = 0.0;    % Hyper diffusivity cutoff ratio
    MU      = 0.0; % Hyperdiffusivity coefficient
    INIT_BLOB = 0; WIPE_TURB = 0; ACT_ON_MODES = 0;
    MU_X    = MU;     %
    MU_Y    = MU;  N_HD    = 4;
    MU_Z    = 2.0; MU_P    = 0.0;     %
    MU_J    = 0.0; LAMBDAD = 0.0;
    NOISE0  = 0.0e-5; % Init noise amplitude
    BCKGD0  = 1.0;    % Init background
GRADB   = 1.0;CURVB   = 1.0;
end
%%-------------------------------------------------------------------------
% RUN
setup
if RUN
    system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 6 ',HELAZDIR,'bin/',EXECNAME,' 2 1 3 0; cd ../../../wk'])
end

% Load results
filename = [SIMID,'/',PARAMS,'/'];
LOCALDIR  = [HELAZDIR,'results/',filename,'/'];
data = compile_results(LOCALDIR,0,0); %Compile the results from first output found to JOBNUMMAX if existing

%linear growth rate (adapted for 2D zpinch and fluxtube)
trange = [0.5 1]*data.Ts3D(end);
nplots = 0;
lg = compute_fluxtube_growth_rate(data,trange,nplots);
[gmax,     kmax] = max(lg.g_ky(:,end));
[gmaxok, kmaxok] = max(lg.g_ky(:,end)./lg.ky);
msg = sprintf('gmax = %2.2f, kmax = %2.2f',gmax,lg.ky(kmax)); disp(msg);
msg = sprintf('gmax/k = %2.2f, kmax/k = %2.2f',gmaxok,lg.ky(kmaxok)); disp(msg);
    
    
    g_max(i) = gmax;
    k_max(i) = kmax;
    
    g_avg(i) = lg.avg_g;
    g_std(i) = lg.std_g;
    
    i = i + 1;
end
%%

% plot(KT_a,max(g_max,0));
y_ = g_avg; 
e_ = g_std;

y_ = y_.*(y_-e_>0);
e_ = e_ .* (y_>0);
errorbar(NU_a,y_,e_,...
    'LineWidth',1.2,...
    'DisplayName',['(',num2str(P),',',num2str(P/2),')']); 
hold on;
title(['$\kappa_T=$',num2str(K_T),' $k_y=$',num2str(ky_),' (CLOS = 0)']);
legend('show'); xlabel('$\nu_{DGGK}$'); ylabel('$\gamma$');
drawnow
end

