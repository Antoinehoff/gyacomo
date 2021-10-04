for CO = [1]
    RUN = 1; % To run or just to load
addpath(genpath('../matlab')) % ... add
default_plots_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
CLUSTER.TIME  = '99:00:00'; % allocation time hh:mm:ss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 0.1;   % Collision frequency
TAU     = 1.0;    % e/i temperature ratio
ETAB    = 1.0;
ETAN    = 1/0.6;    % Density gradient
ETAT    = 0.0;    % Temperature gradient
NU_HYP  = 0.0;   % Hyperdiffusivity coefficient
LAMBDAD = 0.0;
NOISE0  = 1.0e-5; % Init noise amplitude
BCKGD0  = 0.0;    % Init background
%% GRID PARAMETERS
N       = 150;     % Frequency gridpoints (Nkx = N/2)
L       = 100;     % Size of the squared frequency domain
KXEQ0   = 1;      % put kx = 0
MU_P    = 0.0;     % Hermite  hyperdiffusivity -mu_p*(d/dvpar)^4 f
MU_J    = 0.0;     % Laguerre hyperdiffusivity -mu_j*(d/dvperp)^4 f
%% TIME PARMETERS
TMAX    = 400;  % Maximal time unit
DT      = 1e-2;   % Time step
SPS0D   = 1;      % Sampling per time unit for 2D arrays
SPS2D   = 0;      % Sampling per time unit for 2D arrays
SPS3D   = 2;      % Sampling per time unit for 2D arrays
SPS5D   = 2;    % Sampling per time unit for 5D arrays
SPSCP   = 0;    % Sampling per time unit for checkpoints
RESTART = 0;      % To restart from last checkpoint
JOB2LOAD= 00;
%% OPTIONS
% SIMID   = 'v3.6_kobayashi_lin';  % Name of the simulation
% SIMID   = 'v3.2_CO_damping';  % Name of the simulation
% SIMID   = 'CO_Patchwork_damping';  % Name of the simulation
SIMID   = 'test_GF_closure';  % Name of the simulation
% SIMID   = 'v3.2_entropy_mode_linear';  % Name of the simulation
NON_LIN = 0 *(1-KXEQ0);   % activate non-linearity (is cancelled if KXEQ0 = 1)
% Collision operator
% (0 : L.Bernstein, 1 : Dougherty, 2: Sugama, 3 : Pitch angle, 4 : Full Couloumb ; +/- for GK/DK)
% CO      = 1;
INIT_ZF = 0; ZF_AMP = 0.0;
CLOS    = 1;   % Closure model (0: =0 truncation, 1: gyrofluid closure (p+2j<=Pmax))
NL_CLOS = 0;   % nonlinear closure model (0: =0 nmax = jmax, 1: nmax = jmax-j, >1 : nmax = NL_CLOS)
KERN    = 0;   % Kernel model (0 : GK)
INIT_PHI= 0;   % Start simulation with a noisy phi
%% OUTPUTS
W_DOUBLE = 0;
W_GAMMA  = 0;
W_PHI    = 1;
W_NA00   = 1;
W_NAPJ   = 1;
W_SAPJ   = 0;
W_DENS   = 0;
W_TEMP   = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unused
% DK    = 0;  % Drift kinetic model (put every kernel_n to 0 except n=0 to 1)
JOBNUM  = 00;
KPAR    = 0.0;    % Parellel wave vector component
HD_CO   = 0.5;    % Hyper diffusivity cutoff ratio
kmax    = N*pi/L;% Highest fourier mode
MU      = NU_HYP/(HD_CO*kmax)^4; % Hyperdiffusivity coefficient
Nz      = 1;      % number of perpendicular planes (parallel grid)
q0      = 1.0;    % safety factor
shear   = 0.0;    % magnetic shear
eps     = 0.0;    % inverse aspect ratio
INIT_BLOB = 0; WIPE_TURB = 0; WIPE_ZF = 0;
%% PARAMETER SCANS

if 1
%% Parameter scan over PJ
% PA = [2 4];
% JA = [1 2];
PA = [5];
JA = [2];
DTA= DT*ones(size(JA));%./sqrt(JA);
% DTA= DT;
mup_ = MU_P;
muj_ = MU_J;
Nparam = numel(PA);
param_name = 'PJ';
gamma_Ni00 = zeros(Nparam,floor(N/2)+1);
gamma_Nipj = zeros(Nparam,floor(N/2)+1);
gamma_phi  = zeros(Nparam,floor(N/2)+1);
% Ni00_ST  = zeros(Nparam,floor(N/2)+1,TMAX/SPS3D);
%  PHI_ST  = zeros(Nparam,floor(N/2)+1,TMAX/SPS3D);
for i = 1:Nparam
    % Change scan parameter
    PMAXE = PA(i); PMAXI = PA(i);
    JMAXE = JA(i); JMAXI = JA(i);
    DT = DTA(i);
    setup
    % Run linear simulation
    if RUN
        system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 6 ./../../../bin/helaz 1 6; cd ../../../wk'])
%         system(['cd ../results/',SIMID,'/',PARAMS,'/; ./../../../bin/helaz; cd ../../../wk'])
    end
%     Load and process results
    %%
    filename = ['../results/',SIMID,'/',PARAMS,'/outputs_00.h5'];
    load_results
    for ikx = 1:N/2+1
        %find if there is a tail of 0 (highest damped mode)

        tend   = max(Ts3D(abs(Ni00(ikx,1,1,:))~=0));
        tstart   = 0.8*tend;
        [~,itstart] = min(abs(Ts3D-tstart));
        [~,itend]   = min(abs(Ts3D-tend));
        gamma_Ni00(i,ikx) = (LinearFit_s(Ts3D(itstart:itend)',(squeeze(abs(Ni00(ikx,1,1,itstart:itend))))));
        gamma_phi (i,ikx) = (LinearFit_s(Ts3D(itstart:itend)',(squeeze(abs(PHI (ikx,1,1,itstart:itend))))));
%         Ni00_ST(i,ikx,:) = squeeze((Ni00(ikx,1,1,1:TMAX/SPS3D)));
%          PHI_ST(i,ikx,:) = squeeze((PHI (ikx,1,1,1:TMAX/SPS3D)));
    end
    tend   = Ts5D(end); tstart   = 0.4*tend;
    [~,itstart] = min(abs(Ts5D-tstart));
    [~,itend]   = min(abs(Ts5D-tend));
    for ikx = 1:N/2+1
        gamma_Nipj(i,ikx) = LinearFit_s(Ts5D(itstart:itend)',squeeze(max(max(abs(Nipj(:,:,ikx,1,1,itstart:itend)),[],1),[],2)));
    end
    gamma_Ni00(i,:) = real(gamma_Ni00(i,:));% .* (gamma_Ni00(i,:)>=0.0));
    gamma_Nipj(i,:) = real(gamma_Nipj(i,:));% .* (gamma_Nipj(i,:)>=0.0));
    % Clean output
%     system(['rm -r ',BASIC.RESDIR]);
end

if 1
%% Plot
SCALE = 1;%sqrt(2);
fig = figure; FIGNAME = 'linear_study';
plt = @(x) x;
% subplot(211)
    for i = 1:Nparam
        clr       = line_colors(mod(i-1,numel(line_colors(:,1)))+1,:);
        linestyle = line_styles(floor((i-1)/numel(line_colors(:,1)))+1);
        plot(plt(SCALE*kx(1:numel(kx))),plt(gamma_Ni00(i,1:end)),...
            'Color',clr,...
            'LineStyle',linestyle{1},'Marker','^',...
...%             'DisplayName',['$\eta=',num2str(ETAB/ETAN),'$, $\nu_{',CONAME,'}=',num2str(NU),'$, $P=',num2str(PA(i)),'$, $J=',num2str(JA(i)),'$']);
            'DisplayName',[CONAME,', $P,J=',num2str(PA(i)),',',num2str(JA(i)),'$']);
        hold on;
    end
    grid on; xlabel('$k_z\rho_s^{R}$'); ylabel('$\gamma(N_i^{00})L_\perp/c_s$'); xlim([0.0,max(kx)]);
    title(['$\eta=',num2str(ETAB/ETAN),'$, $\nu_{',CONAME,'}=',num2str(NU),'$'])
%     title(['$\nabla N = 0$', ', $\nu=',num2str(NU),'$'])
    legend('show'); %xlim([0.01,10])
saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.fig']);
saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.png']);
end
end
if 0
%% Space time
    [YT,XT] = meshgrid(Ts3D,kx);
    figure;
%     pclr = surf(XT,YT,squeeze(abs(PHI_ST(1,:,:)))); set(pclr, 'edgecolor','none'); colorbar;
%     pclr = pcolor(XT,YT,squeeze(abs(Ni00_ST(1,:,:)))); set(pclr, 'edgecolor','none'); colorbar;
    semilogy(Ts3D(1:TMAX/SPS3D),squeeze(abs(PHI_ST(1,50:5:100,:))));
end
if 0 
%% Trajectories of some modes
figure;
for i = 1:10:N/2+1
    semilogy(Ts3D,squeeze(abs(Ne00(i,1,1,:))),'DisplayName',['k=',num2str(kx(i))]); hold on;
end
end
end
