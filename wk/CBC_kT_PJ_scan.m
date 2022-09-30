addpath(genpath('../matlab')) % ... add
default_plots_options
HELAZDIR = '/home/ahoffman/HeLaZ/';
EXECNAME = 'helaz3';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KT_a = [3:0.5:5];
P_a_6   = [6 6 6 6 6 6 6];
J_a_6   = [0 1 2 3 4 5 6];
P_a_10  = 10*ones(1,6);
J_a_10  = 5:10;
% P_a_10b = 10*ones(1,7);
% J_a_10b = 10:17;
P_a_20  = 20*ones(1,11);
J_a_20  = 10:20;

P_a = [P_a_20]; J_a = [J_a_20];
% KT_a = 5.0; P_a = 20; J_a = 20;

g_max= zeros(numel(P_a),numel(KT_a));
g_avg= g_max*0;
g_std= g_max*0;
k_max= g_max*0;

CO    = 'DG'; GKCO = 0;
NU    = 0.0;
DT    = 7e-3;
TMAX  = 40;
ky_   = 0.15;
SIMID = 'linear_CBC_kT_threshold';  % Name of the simulation
RUN   = 0;


for i = 1:numel(P_a)
P = P_a(i); J = J_a(i);
    j=1;
    %Set Up parameters
    for K_Ti = KT_a
        CLUSTER.TIME  = '99:00:00'; % allocation time hh:mm:ss
        TAU     = 1.0;    % e/i temperature ratio
        K_Ni = 2.22; K_Ne = K_Ni;
        K_Te     = K_Ti;            % Temperature '''
        SIGMA_E = 0.0233380;   % mass ratio sqrt(m_a/m_i) (correct = 0.0233380)
        KIN_E   = 0;     % 1: kinetic electrons, 2: adiabatic electrons
        BETA    = 0e-1;     % electron plasma beta 
        PMAXE   = P; JMAXE   = J;
        PMAXI   = P; JMAXI   = J;
        NX      = 8;    % real space x-gridpoints
        NY      = 6;     %     ''     y-gridpoints
        LX      = 2*pi/0.15;   % Size of the squared frequency domain
        LY      = 2*pi/ky_;
        NZ      = 24;    % number of perpendicular planes (parallel grid)
        NPOL    = 1; SG = 0; NEXC = 1;
        GEOMETRY= 's-alpha';
        Q0      = 1.4;    % safety factor
        SHEAR   = 0.8;    % magnetic shear (Not implemented yet)
        EPS     = 0.18;    % inverse aspect ratio
        SPS0D = 1; SPS2D = 0; SPS3D = 5;SPS5D= 1/5; SPSCP = 0;
        JOB2LOAD= -1;
        LINEARITY = 'linear';   % activate non-linearity (is cancelled if KXEQ0 = 1)
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
        MU_Z    = 1.0; MU_P    = 0.0;     %
        MU_J    = 0.0; LAMBDAD = 0.0;
        NOISE0  = 1.0e-4; % Init noise amplitude
        BCKGD0  = 0.0;    % Init background
        GRADB   = 1.0;CURVB   = 1.0;

        %%-------------------------------------------------------------------------
        % RUN
        setup
        if RUN
            system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 4 ',HELAZDIR,'bin/',EXECNAME,' 2 2 1 0; cd ../../../wk'])
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

        g_max(i,j) = gmax;
        k_max(i,j) = kmax;
        
        [g_avg(i,j), ik_] = max(lg.avg_g);
        g_std(i,j) = max(lg.std_g(ik_));

        j = j + 1;
        if 0
        %% Verify gamma time trace
        figure
        for ik_ = 1:numel(lg.ky)
            plot(lg.trange(2:end),lg.g_ky(ik_,:)','DisplayName',['$k_y=',num2str(lg.ky(ik_)),'$']); hold on;
        end
        xlabel('$t$'); ylabel('$\gamma$');
        title(data.param_title); legend('show');
        drawnow
        end
    end
end

if 1
%% PLOTS
ERR_WEIGHT = 1/3; %weight of the error to compute marginal stability
%% Superposed 1D plots
sz_ = size(g_max);
figure
for i = 1:sz_(1)
    y_ = g_avg(i,:); 
    e_ = g_std(i,:);

    y_ = y_.*(y_-e_*ERR_WEIGHT>0);
    e_ = e_ .* (y_>0);
        errorbar(KT_a,y_,e_,...
            'LineWidth',1.2,...
            'DisplayName',['(',num2str(P_a(i)),',',num2str(J_a(i)),')']); 
%     plot(KT_a,y_,...
%         'LineWidth',1.2,...
%         'DisplayName',['(',num2str(P_a(i)),',',num2str(J_a(i)),')']); 
    hold on;

end
title('Linear CBC $K_T$ threshold');
legend('show'); xlabel('$K_T$'); ylabel('$\max_{k_y}(\gamma_k)$');
drawnow
%% Color map
[NP__, KT__] = meshgrid(P_a+2*J_a, KT_a);
% GG_ = g_avg;
GG_ = g_avg .* (g_avg-g_std*ERR_WEIGHT > 0);
figure;
% pclr = pcolor(KT__,NP__,g_max'); set(pclr,'EdgeColor','none');
pclr = imagesc(KT_a,1:numel(P_a),GG_);
LABELS = [];
for i_ = 1:numel(P_a)
    LABELS = [LABELS; '(',sprintf('%2.0f',P_a(i_)),',',sprintf('%2.0f',J_a(i_)),')'];
end
yticks(1:numel(P_a));
yticklabels(LABELS);
xlabel('$\kappa_T$'); ylabel('$(P,J)$');
title('Linear ITG threshold in CBC');
colormap(bluewhitered);
%%
%%
end

