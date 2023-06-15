gyacomodir  = pwd;
gyacomodir = gyacomodir(1:end-2);
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add% EXECNAME = 'gyacomo_1.0';
% EXECNAME = 'gyacomo23_sp';
EXECNAME = 'gyacomo23_dp';
% EXECNAME = 'gyacomo23_debug';
CLUSTER.TIME  = '99:00:00'; % allocation time hh:mm:ss
%%
SIMID = 'p2_linear_new';  % Name of the simulation
RERUN   = 0; % rerun if the  data does not exist
RUN     = 0;
K_T = 5.3;
P_a = [2 4 8 16 32 64];
% P_a = 32;
ky_a= [0.05:0.05:1.0];
% ky_a= [0.25:0.05:0.65];
% collision setting
CO        = 'DG';
GKCO      = 0; % gyrokinetic operator
COLL_KCUT = 1.75;
NU        = 1e-3;
% model
KIN_E   = 0;         % 1: kinetic electrons, 2: adiabatic electrons
BETA    = 1e-4;     % electron plasma beta
% background gradients setting
K_N    = 2.22;            % Density '''
% Geometry
% GEOMETRY= 'miller';
GEOMETRY= 's-alpha';
SHEAR   = 0.8;    % magnetic shear
% time and numerical grid
DT0    = 1e-3;
TMAX   = 30;
% arrays for the result
g_ky = zeros(numel(ky_a),numel(P_a),1);
g_avg= g_ky*0;
g_std= g_ky*0;
% Naming of the collision operator
if GKCO
    CONAME = [CO,'GK'];
else
    CONAME = [CO,'DK'];
end
j = 1;
for P = P_a
J = P/2;
i = 1;
for ky = ky_a
    %% PHYSICAL PARAMETERS
    TAU     = 1.0;            % e/i temperature ratio
    % SIGMA_E = 0.05196152422706632;   % mass ratio sqrt(m_a/m_i) (correct = 0.0233380)
    SIGMA_E = 0.0233380;   % mass ratio sqrt(m_a/m_i) (correct = 0.0233380)
    K_Te    = K_T;            % ele Temperature '''
    K_Ti    = K_T;            % ion Temperature '''
    K_Ne    = K_N;            % ele Density '''
    K_Ni    = K_N;            % ion Density gradient drive
    %% GRID PARAMETERS
    DT = DT0/sqrt(J);
    PMAX   = P;     % Hermite basis size
    JMAX   = P/2;     % Laguerre "
    NX      = 8;    % real space x-gridpoints
    NY      = 2;
    LX      = 2*pi/0.8;   % Size of the squared frequency domain
    LY      = 2*pi/ky;     % Size of the squared frequency domain
    NZ      = 24;    % number of perpendicular planes (parallel grid)
    NPOL    = 1;
    SG      = 0;     % Staggered z grids option
    NEXC    = 1;     % To extend Lx if needed (Lx = Nexc/(kymin*shear))
    %% GEOMETRY
    % GEOMETRY= 's-alpha';
    EPS     = 0.18;   % inverse aspect ratio
    Q0      = 1.4;    % safety factor
    KAPPA   = 1.0;    % elongation
    DELTA   = 0.0;    % triangularity
    ZETA    = 0.0;    % squareness
    PARALLEL_BC = 'dirichlet'; %'dirichlet','periodic','shearless','disconnected'
%     PARALLEL_BC = 'periodic'; %'dirichlet','periodic','shearless','disconnected'
    SHIFT_Y = 0.0;
    %% TIME PARMETERS
    DTSAVE0D = 1;      % Sampling per time unit for 0D arrays
    DTSAVE2D = -1;     % Sampling per time unit for 2D arrays
    DTSAVE3D = 2;      % Sampling per time unit for 3D arrays
    DTSAVE5D = 100;     % Sampling per time unit for 5D arrays
    JOB2LOAD = -1;     % Start a new simulation serie
    %% OPTIONS
    LINEARITY = 'linear';   % activate non-linearity (is cancelled if KXEQ0 = 1)
    % Collision operator
    ABCO    = 1; % INTERSPECIES collisions
    INIT_ZF = 0; ZF_AMP = 0.0;
    CLOS    = 0;   % Closure model (0: =0 truncation, 1: v^Nmax closure (p+2j<=Pmax))s
    NL_CLOS = 0;   % nonlinear closure model (-2:nmax=jmax; -1:nmax=jmax-j; >=0:nmax=NL_CLOS)
    KERN    = 0;   % Kernel model (0 : GK)
    INIT_OPT= 'phi';   % Start simulation with a noisy mom00/phi/allmom
    NUMERICAL_SCHEME = 'RK4'; % RK2,SSPx_RK2,RK3,SSP_RK3,SSPx_RK3,IMEX_SSP2,ARK2,RK4,DOPRI5
    %% OUTPUTS
    W_DOUBLE = 0;
    W_GAMMA  = 1; W_HF     = 1;
    W_PHI    = 1; W_NA00   = 1;
    W_DENS   = 0; W_TEMP   = 1;
    W_NAPJ   = 0; W_SAPJ   = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % unused
    NA      = 1;
    ADIAB_E = (NA==1);
    HD_CO   = 0.0;    % Hyper diffusivity cutoff ratio
    MU      = 0.0; % Hyperdiffusivity coefficient
    INIT_BLOB = 0; WIPE_TURB = 0; ACT_ON_MODES = 0;
    MU_X    = MU;     %
    MU_Y    = MU;     %
    N_HD    = 4;
    HYP_V   = 'none';
    HRCY_CLOS = 'truncation';   % Closure model for higher order moments
    DMAX      = -1;
    NLIN_CLOS = 'truncation';   % Nonlinear closure model for higher order moments
    NMAX      = 0;
    MU_Z    = 1.0;     %
    MU_P    = 0.0;     %
    MU_J    = 0.0;     %
    LAMBDAD = 0.0;
    NOISE0  = 1.0e-5; % Init noise amplitude
    BCKGD0  = 0.0;    % Init background
    k_gB   = 1.0;
    k_cB   = 1.0;
    %% RUN
    setup
    % naming
    filename = [SIMID,'/',PARAMS,'/'];
    LOCALDIR  = [gyacomodir,'results/',filename,'/'];
    % check if data exist to run if no data
    data_ = {};
    try
        data_ = compile_results_low_mem(data_,LOCALDIR,00,00);
        Ntime = numel(data_.Ts0D);
    catch
        data_.outfilenames = [];
    end
    if RUN && (RERUN || isempty(data_.outfilenames) || Ntime < 10)
        % system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 2 ',gyacomodir,'bin/',EXECNAME,' 1 2 1 0; cd ../../../wk'])
        system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 4 ',gyacomodir,'bin/',EXECNAME,' 2 2 1 0; cd ../../../wk'])
        % system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 6 ',gyacomodir,'bin/',EXECNAME,' 3 2 1 0; cd ../../../wk'])
    end
    data_    = compile_results_low_mem(data_,LOCALDIR,00,00);
    [data_.PHI, data_.Ts3D] = compile_results_3D(LOCALDIR,00,00,'phi');
    if numel(data_.Ts3D)>10
        if numel(data_.Ts3D)>10
        % Load results after trying to run
        filename = [SIMID,'/',PARAMS,'/'];
        LOCALDIR  = [gyacomodir,'results/',filename,'/'];

        data_    = compile_results_low_mem(data_,LOCALDIR,00,00);
        [data_.PHI, data_.Ts3D] = compile_results_3D(LOCALDIR,00,00,'phi');

        % linear growth rate (adapted for 2D zpinch and fluxtube)
        options.TRANGE = [0.5 1]*data_.Ts3D(end);
        options.NPLOTS = 0; % 1 for only growth rate and error, 2 for omega local evolution, 3 for plot according to z
        options.GOK    = 0; %plot 0: gamma 1: gamma/k 2: gamma^2/k^3

        [~,it1] = min(abs(data_.Ts3D-0.5*data_.Ts3D(end))); % start of the measurement time window
        [~,it2] = min(abs(data_.Ts3D-1.0*data_.Ts3D(end))); % end of ...
        field   = 0;
        field_t = 0;
        for ik = 2:NY/2+1
            field   = squeeze(sum(abs(data_.PHI),3)); % take the sum over z
            field_t = squeeze(field(ik,1,:)); % take the kx =0, ky = ky mode only
            to_measure  = log(field_t(it1:it2));
            tw = double(data_.Ts3D(it1:it2));
    %         gr = polyfit(tw,to_measure,1);
            gr = fit(tw,to_measure,'poly1');
            err= confint(gr);
            g_ky(i,j,ik)  = gr.p1;
            g_std(i,j,ik) = abs(err(2,1)-err(1,1))/2;
        end
        [gmax, ikmax] = max(g_ky(i,j,:));

        msg = sprintf('gmax = %2.2f, kmax = %2.2f',gmax,data_.grids.ky(ikmax)); disp(msg);
        end
    end
    
    i = i + 1;
end
j = j + 1;
end

if 0
%% Check time evolution
figure;
to_measure  = log(field_t);
plot(data_.Ts3D,to_measure); hold on
plot(data_.Ts3D(it1:it2),to_measure(it1:it2),'--');
end

%% take max growth rate among z coordinate
y_ = g_ky(:,:,2); 
e_ = g_std(:,:,2);

%%
if(numel(ky_a)>1 && numel(P_a)>1)
%% Save metadata
 pmin  = num2str(min(P_a));   pmax = num2str(max(P_a));
 kymin = num2str(min(ky_a));  kymax= num2str(max(ky_a));
filename = [num2str(NX),'x',num2str(NZ),'_ky_',kymin,'_',kymax,...
            '_P_',pmin,'_',pmax,'_',CONAME,'_',num2str(NU),'_kT_',num2str(K_T),'.mat'];
metadata.name   = filename;
metadata.kymin  = ky;
metadata.title  = ['$\nu_{',CONAME,'}=$',num2str(NU),'$\kappa_T=$',num2str(K_T),', $\kappa_N=$',num2str(K_N)];
metadata.par    = [num2str(NX),'x1x',num2str(NZ)];
metadata.nscan  = 2;
metadata.s2name = '$P$';
metadata.s2     = P_a;
metadata.s1name = '$ky$';
metadata.s1     = ky_a;
metadata.dname  = '$\gamma c_s/R$';
metadata.data   = y_;
metadata.err    = e_;
save([SIMDIR,filename],'-struct','metadata');
disp(['saved in ',SIMDIR,filename]);
clear metadata tosave
end
