gyacomodir = '/home/ahoffman/gyacomo/';
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add% EXECNAME = 'gyacomo_1.0';
% EXECNAME = 'gyacomo_dbg';
EXECNAME = 'gyacomo';
%%
NU_a = [0:0.01:0.05]*0.1/0.045;
% P_a  = [2 4 6 8 10 12 16];
% NU_a = 0.1;
P_a  = [2 4];


CO      = 'DG';
COLL_KCUT = 1.75;

K_Ti  = 8.0;
K_Ni  = 2.0;
DT    = 1e-3;
TMAX  = 20;
kymin = 0.35;
Ny    = 2;
SIMID = 'Ajay_CH4_lin_ITG';  % Name of the simulation
RUN   = 1;

g_ky = zeros(numel(NU_a),numel(P_a),Ny/2);
g_avg= g_ky*0;
g_std= g_ky*0;

j = 1;
for P = P_a

i = 1;
for NU = NU_a
    %Set Up parameters
    CLUSTER.TIME  = '99:00:00'; % allocation time hh:mm:ss
    TAU     = 1.0;    % e/i temperature ratio
    K_Ni = 2.0; K_Ne = K_Ni;
    K_Te     = K_Ti;            % Temperature '''
    SIGMA_E = 0.0233380;   % mass ratio sqrt(m_a/m_i) (correct = 0.0233380)
    KIN_E   = 0;     % 1: kinetic electrons, 2: adiabatic electrons
    BETA    = 1e-3;     % electron plasma beta 
    J = P/2;
%     J = 2;
    PMAXE   = P; JMAXE   = J;
    PMAXI   = P; JMAXI   = J;
    NX      = 16;    % real space x-gridpoints
    NY      = Ny;     %     ''     y-gridpoints
    LX      = 142.9;   % Size of the squared frequency domain
    LY      = 2*pi/kymin;
    NZ      = 24;    % number of perpendicular planes (parallel grid)
    NPOL    = 1; SG = 0;
%     GEOMETRY= 's-alpha';
    GEOMETRY= 'miller';
    EPS     = 0.18;   % inverse aspect ratio
    Q0      = 1.4;    % safety factor
    SHEAR   = 0.8;    % magnetic shear
    KAPPA   = 1.0;    % elongation
    DELTA   = 0.0;    % triangularity
    ZETA    = 0.0;    % squareness
    NEXC    = 1;      % To extend Lx if needed (Lx = Nexc/(kymin*shear))
    SPS0D = 1; SPS2D = -1; SPS3D = 1;SPS5D= 1/5; SPSCP = 0;
    JOB2LOAD= -1;
    LINEARITY = 'linear';   % activate non-linearity (is cancelled if KXEQ0 = 1)
    GKCO    = 1; % gyrokinetic operator
    ABCO    = 1; % interspecies collisions
    INIT_ZF = 0; ZF_AMP = 0.0;
    CLOS    = 0;   % Closure model (0: =0 truncation, 1: v^Nmax closure (p+2j<=Pmax))s
    NL_CLOS = 0;   % nonlinear closure model (-2:nmax=jmax; -1:nmax=jmax-j; >=0:nmax=NL_CLOS)
    KERN    = 0;   % Kernel model (0 : GK)
    INIT_OPT= 'mom00';   % Start simulation with a noisy mom00/phi/allmom
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
    NOISE0  = 1.0e-5; % Init noise amplitude
    BCKGD0  = 0.0;    % Init background
    GRADB   = 1.0;CURVB   = 1.0;
    %%-------------------------------------------------------------------------
    % RUN
    setup
    if RUN
%         system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 6 ',gyacomodir,'bin/',EXECNAME,' 1 6 1 0; cd ../../../wk'])
        system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 4 ',gyacomodir,'bin/',EXECNAME,' 1 1 4 0; cd ../../../wk'])
    end

    % Load results
    filename = [SIMID,'/',PARAMS,'/'];
    LOCALDIR  = [gyacomodir,'results/',filename,'/'];
    data = compile_results(LOCALDIR,0,0); %Compile the results from first output found to JOBNUMMAX if existing

    % linear growth rate (adapted for 2D zpinch and fluxtube)
    options.TRANGE = [0.5 1]*data.Ts3D(end);
    options.NPLOTS = 0; % 1 for only growth rate and error, 2 for omega local evolution, 3 for plot according to z
    options.GOK    = 0; %plot 0: gamma 1: gamma/k 2: gamma^2/k^3
    lg = compute_fluxtube_growth_rate(data,options);
    [gmax,     kmax] = max(lg.g_ky(:,end));
    [gmaxok, kmaxok] = max(lg.g_ky(:,end)./lg.ky);
    msg = sprintf('gmax = %2.2f, kmax = %2.2f',gmax,lg.ky(kmax)); disp(msg);
    msg = sprintf('gmax/k = %2.2f, kmax/k = %2.2f',gmaxok,lg.ky(kmaxok)); disp(msg);

    
    g_ky(i,j,:)  = lg.avg_g;
    
    g_avg(i,j,:) = lg.avg_g;
    g_std(i,j,:) = lg.std_g;
    
    i = i + 1;
end
j = j + 1;
end

if 1
%% Study of the peak growth rate
figure

y_ = g_avg; 
e_ = g_std;

% filter to noisy data
y_ = y_.*(y_-e_>0);
e_ = e_ .* (y_>0);

[y_,idx_] = max(g_avg,[],3); 
for i = 1:numel(idx_)
    e_ = g_std(:,:,idx_(i));
end

colors_ = summer(numel(NU_a));
subplot(121)
for i = 1:numel(NU_a)
    errorbar(P_a,y_(i,:),e_(i,:),...
        'LineWidth',1.2,...
        'DisplayName',['$\nu=$',num2str(NU_a(i))],...
        'color',colors_(i,:)); 
    hold on;
end
title(['$\kappa_T=$',num2str(K_Ti),' $k_y=k_y^{max}$']);
legend('show'); xlabel('$P$, $J=P/2$'); ylabel('$\gamma$');

colors_ = jet(numel(P_a));
subplot(122)
for j = 1:numel(P_a)
    errorbar(NU_a,y_(:,j),e_(:,j),...
        'LineWidth',1.2,...
        'DisplayName',['(',num2str(P_a(j)),',',num2str(P_a(j)/2),')'],...
        'color',colors_(j,:)); 
    hold on;
end
title(['$\kappa_T=$',num2str(K_Ti),' $k_y=k_y^{max}$']);
legend('show'); xlabel(['$\nu_{',CO,'}$']); ylabel('$\gamma$');
end

if 0
%% Pcolor of the peak
figure;
[XX_,YY_] = meshgrid(NU_a,P_a);
pclr=pcolor(XX_,YY_,y_'); set(pclr,'EdgeColor','none');
end