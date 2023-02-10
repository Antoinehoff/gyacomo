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
SIMID = 'p2_CBC_convergence_KT_PJ';  % Name of the simulation
% SIMID = 'dbg';  % Name of the simulation
RERUN   = 0; % rerun if the data does not exist
RUN     = 0;
KT_a = [3:0.5:7];
% KT_a = [3];
% P_a  = [22];
P_a  = [2:2:30];
% P_a  = [2:12];
J_a  = floor(P_a/2);
% collision setting
CO        = 'DG';
NU        = 0.05;
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
DT    = 1e-3;
TMAX  = 25;
kymin = 0.3;
NY    = 2;
% arrays for the result
g_ky = zeros(numel(KT_a),numel(P_a),NY/2+1);
g_avg= g_ky*0;
g_std= g_ky*0;

j = 1;
for P = P_a
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
%     P = 20;
    J = floor(P/2);
    PMAXE   = P;     % Hermite basis size of electrons
    JMAXE   = J;     % Laguerre "
    PMAXI   = P;     % " ions
    JMAXI   = J;     % "
    NX      = 12;    % real space x-gridpoints
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
    ABCO    = 1; % interspecies collisions
    INIT_ZF = 0; ZF_AMP = 0.0;
    CLOS    = 0;   % Closure model (0: =0 truncation, 1: v^Nmax closure (p+2j<=Pmax))s
    NL_CLOS = 0;   % nonlinear closure model (-2:nmax=jmax; -1:nmax=jmax-j; >=0:nmax=NL_CLOS)
    KERN    = 0;   % Kernel model (0 : GK)
    INIT_OPT= 'phi';   % Start simulation with a noisy mom00/phi/allmom
    NUMERICAL_SCHEME = 'RK4'; % RK2,SSPx_RK2,RK3,SSP_RK3,SSPx_RK3,IMEX_SSP2,ARK2,RK4,DOPRI5
    %% OUTPUTS
    W_DOUBLE = 1;
    W_GAMMA  = 1; W_HF     = 1;
    W_PHI    = 1; W_NA00   = 1;
    W_DENS   = 1; W_TEMP   = 1;
    W_NAPJ   = 1; W_SAPJ   = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % unused
    HD_CO   = 0.0;    % Hyper diffusivity cutoff ratio
    MU      = 0.0; % Hyperdiffusivity coefficient
    INIT_BLOB = 0; WIPE_TURB = 0; ACT_ON_MODES = 0;
    MU_X    = MU;     %
    MU_Y    = MU;     %
    N_HD    = 4;
    MU_Z    = 0.2;     %
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

    
    i = i + 1;
end
j = j + 1;
end

if 0 
%% Check time evolution
figure;
plot(data_.Ts3D,to_measure); hold on
plot(data_.Ts3D(it1:it2),to_measure(it1:it2),'--');
end

if 1
%% Study of the peak growth rate
figure

y_ = g_ky; 
e_ = 0.05;

% filter growth rate with less than 0.05 value
y_ = y_.*(y_-e_>0);
e_ = e_ .* (y_>0);

[y_,idx_] = max(g_ky,[],3); 
for i = 1:numel(idx_)
    e_ = g_std(:,:,idx_(i));
end

colors_ = jet(numel(KT_a));
subplot(121)
for i = 1:numel(KT_a)
%     errorbar(P_a,y_(i,:),e_(i,:),...
%         'LineWidth',1.2,...
%         'DisplayName',['$\nu=$',num2str(KT_a(i))],...
%         'color',colors_(i,:)); 
    plot(P_a,y_(i,:),'s-',...
        'LineWidth',2.0,...
        'DisplayName',['$\kappa_T=$',num2str(KT_a(i))],...
        'color',colors_(i,:)); 
    hold on;
end
title(['$\nu_{',CO,'}=$',num2str(NU),', $\kappa_N=$',num2str(K_N),', $k_y=$',num2str(kymin)]);
legend('show'); xlabel('$P$, $J=P/2$'); ylabel('$\gamma$');

colors_ = jet(numel(P_a));
subplot(122)
for j = 1:numel(P_a)
% errorbar(KT_a,y_(:,j),e_(:,j),...
%     'LineWidth',1.2,...
%     'DisplayName',['(',num2str(P_a(j)),',',num2str(J_a(j)),')'],...
%     'color',colors_(j,:)); 
    plot(KT_a,y_(:,j),'s-',...
        'LineWidth',2.0,...
        'DisplayName',['(',num2str(P_a(j)),',',num2str(J_a(j)),')'],...
        'color',colors_(j,:)); 
    hold on;
end
title(['$\nu_{',CO,'}=$',num2str(NU),', $\kappa_N=$',num2str(K_N),', $k_y=$',num2str(kymin)]);
legend('show'); xlabel('$\kappa_T$'); ylabel('$\gamma$');
end

if 0
%% Pcolor of the peak
figure;
[XX_,YY_] = meshgrid(KT_a,P_a);
% pclr=pcolor(XX_,YY_,y_'); set(pclr,'EdgeColor','none'); axis ij;
pclr=imagesc_custom(XX_,YY_,y_'.*(y_>0)');
title(['$\nu_{',data.CO,'}=$',num2str(data.NU),', $\kappa_N=$',num2str(K_N),', $k_y=$',num2str(kymin)]);
xlabel('$\kappa_T$'); ylabel('$P$, $J=P/2$');
colormap(bluewhitered)
clb=colorbar; 
clb.Label.String = '$\gamma c_s/R$';
clb.Label.Interpreter = 'latex';
clb.Label.FontSize= 18;
end
%%
%% Save metadata
ktmin = num2str(min(KT_a)); ktmax = num2str(max(KT_a));
 pmin = num2str(min(P_a));   pmax = num2str(max(P_a));
filename = [num2str(NX),'x',num2str(NZ),'_ky_',num2str(kymin),...
            '_kT_',ktmin,'_',ktmax,'_',...
            '_P_',pmin,'_',pmax,'_',data_.CO,'_',num2str(NU),'.mat'];
metadata.name   = filename;
metadata.kymin  = kymin;
metadata.title  = ['$\nu_{',data_.CO,'}=$',num2str(NU),', $\kappa_N=$',num2str(K_N),', $k_y=$',num2str(kymin)];
metadata.par    = data_.PARAMS;
metadata.nscan  = 2;
metadata.s1name = '$\kappa_T$';
metadata.s1     = KT_a;
metadata.s2name = '$P$, $J=P/2$';
metadata.s2     = P_a;
metadata.dname  = '$\gamma c_s/R$';
metadata.data   = y_;
metadata.err    = e_;
metadata.input_file = h5read([data_.localdir,'/outputs_00.h5'],'/files/STDIN.00');
metadata.date   = date;
% tosave.data     = metadata;
save([SIMDIR,filename],'-struct','metadata');
disp(['saved in ',SIMDIR,filename]);
clear metadata tosave