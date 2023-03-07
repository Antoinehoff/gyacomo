gyacomodir  = pwd;
gyacomodir = gyacomodir(1:end-2);
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add% EXECNAME = 'gyacomo_1.0';
% EXECNAME = 'gyacomo_debug';
EXECNAME = 'gyacomo';
CLUSTER.TIME  = '99:00:00'; % allocation time hh:mm:ss
%%
SIMID = 'p2_linear';  % Name of the simulation
RERUN   = 0; % rerun if the data does not exist
RUN     = 1;
KT_a = [3:0.5:6.5 6.96];
NU_a = [0 0.01:0.01:0.1 0.2 0.5 1.0];
% KT_a = 3.5;
% NU_a = 0.5;
P    = 16;
J    = 8;
% collision setting
CO        = 'DG';
GKCO      = 0; % gyrokinetic operator
COLL_KCUT = 1.75;
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
DT     = 2e-3;
TMAX   = 30;
kymin  = 0.3;
NY     = 2;
% arrays for the result
g_ky = zeros(numel(KT_a),numel(NU_a),NY/2+1);
g_avg= g_ky*0;
g_std= g_ky*0;
% Naming of the collision operator
if GKCO
    CONAME = [CO,'GK'];
else
    CONAME = [CO,'DK'];
end
    
j = 1;
for NU = NU_a
i = 1;
for KT = KT_a
    %% PHYSICAL PARAMETERS
    TAU     = 1.0;            % e/i temperature ratio
    % SIGMA_E = 0.05196152422706632;   % mass ratio sqrt(m_a/m_i) (correct = 0.0233380)
    SIGMA_E = 0.0233380;   % mass ratio sqrt(m_a/m_i) (correct = 0.0233380)
    K_Te    = KT;            % ele Temperature '''
    K_Ti    = KT;            % ion Temperature '''
    K_Ne    = K_N;            % ele Density '''
    K_Ni    = K_N;            % ion Density gradient drive
    %% GRID PARAMETERS
    PMAXE   = P;     % Hermite basis size of electrons
    JMAXE   = J;     % Laguerre "
    PMAXI   = P;     % " ions
    JMAXI   = J;     % "
    NX      = 8;    % real space x-gridpoints
    LX      = 2*pi/0.8;   % Size of the squared frequency domain
    LY      = 2*pi/kymin;     % Size of the squared frequency domain
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
    SPS0D   = 1;      % Sampling per time unit for 2D arrays
    SPS2D   = -1;      % Sampling per time unit for 2D arrays
    SPS3D   = 1;      % Sampling per time unit for 2D arrays
    SPS5D   = 1/2;    % Sampling per time unit for 5D arrays
    SPSCP   = 0;    % Sampling per time unit for checkpoints
    JOB2LOAD= -1;
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
    HD_CO   = 0.0;    % Hyper diffusivity cutoff ratio
    MU      = 0.0; % Hyperdiffusivity coefficient
    INIT_BLOB = 0; WIPE_TURB = 0; ACT_ON_MODES = 0;
    MU_X    = MU;     %
    MU_Y    = MU;     %
    N_HD    = 4;
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
    data_ = compile_results(LOCALDIR,0,0); %Compile the results from first output found to JOBNUMMAX if existing
    if RUN && (RERUN || isempty(data_.NU_EVOL) || numel(data_.Ts3D)<10)
        system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 4 ',gyacomodir,'bin/',EXECNAME,' 1 2 2 0; cd ../../../wk'])
%         system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 6 ',gyacomodir,'bin/',EXECNAME,' 3 2 1 0; cd ../../../wk'])
    end
    if ~isempty(data_.NU_EVOL)
        if numel(data_.Ts3D)>10
        % Load results after trying to run
        filename = [SIMID,'/',PARAMS,'/'];
        LOCALDIR  = [gyacomodir,'results/',filename,'/'];

        data_ = compile_results(LOCALDIR,0,0); %Compile the results from first output found to JOBNUMMAX if existing

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
            tw = data_.Ts3D(it1:it2);
    %         gr = polyfit(tw,to_measure,1);
            gr = fit(tw,to_measure,'poly1');
            err= confint(gr);
            g_ky(i,j,ik)  = gr.p1;
            g_std(i,j,ik) = abs(err(2,1)-err(1,1))/2;
        end
        [gmax, ikmax] = max(g_ky(i,j,:));

        msg = sprintf('gmax = %2.2f, kmax = %2.2f',gmax,data_.ky(ikmax)); disp(msg);
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

if(numel(KT_a)>1 && numel(NU_a)>1)
%% Save metadata
ktmin = num2str(min(KT_a)); ktmax = num2str(max(KT_a));
numin = num2str(min(NU_a)); numax = num2str(max(NU_a));
filename = [num2str(NX),'x',num2str(NZ),'_ky_',num2str(kymin),...
            '_P_',num2str(P),'_J_',num2str(J),...
            '_kT_',ktmin,'_',ktmax,...
            '_nu_',numin,'_',numax,'_',CONAME,'.mat'];
metadata.name   = filename;
metadata.kymin  = kymin;
metadata.title  = ['P=',num2str(P),' J=',num2str(J),'$\kappa_N=$',num2str(K_N),', $k_y=$',num2str(kymin)];
metadata.par    = [num2str(NX),'x1x',num2str(NZ)];
metadata.nscan  = 2;
metadata.s1name = '$\kappa_T$';
metadata.s1     = KT_a;
metadata.s2name = ['$\nu_{',CONAME,'}$'];
metadata.s2     = NU_a;
metadata.dname  = '$\gamma c_s/R$';
metadata.data   = y_;
metadata.err    = e_;
metadata.date   = date;
save([SIMDIR,filename],'-struct','metadata');
disp(['saved in ',SIMDIR,filename]);
clear metadata tosave
end