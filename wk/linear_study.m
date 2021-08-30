%clear all;
addpath(genpath('../matlab')) % ... add
default_plots_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
CLUSTER.TIME  = '99:00:00'; % allocation time hh:mm:ss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
% NU      = 1.0;   % Collision frequency
TAU     = 1.0;    % e/i temperature ratio
ETAB    = 1.0;
% ETAN    = ETAN;    % Density gradient
ETAT    = 0.0;    % Temperature gradient
NU_HYP  = 0.0;   % Hyperdiffusivity coefficient
LAMBDAD = 0.0;
NOISE0  = 1.0e-5;
%% GRID PARAMETERS
N       = 50;     % Frequency gridpoints (Nkx = N/2)
L       = 75;     % Size of the squared frequency domain
KXEQ0   = 1;      % put kx = 0
MU_P    = 0.0;     % Hermite  hyperdiffusivity -mu_p*(d/dvpar)^4 f
MU_J    = 0.0;     % Laguerre hyperdiffusivity -mu_j*(d/dvperp)^4 f
%% TIME PARMETERS
TMAX    = 200;  % Maximal time unit
DT      = 1e-2;   % Time step
SPS0D   = 1;      % Sampling per time unit for 2D arrays
SPS2D   = 1;      % Sampling per time unit for 2D arrays
SPS5D   = 1/5;    % Sampling per time unit for 5D arrays
SPSCP   = 0;    % Sampling per time unit for checkpoints
RESTART = 0;      % To restart from last checkpoint
JOB2LOAD= 00;
%% OPTIONS
SIMID   = 'v3.1_lin_analysis';  % Name of the simulation
NON_LIN = 0 *(1-KXEQ0);   % activate non-linearity (is cancelled if KXEQ0 = 1)
% Collision operator
% (0 : L.Bernstein, 1 : Dougherty, 2: Sugama, 3 : Pitch angle, 4 : Full Couloumb ; +/- for GK/DK)
% CO      = 2;
INIT_ZF = 0; ZF_AMP = 0.0;
CLOS    = 0;   % Closure model (0: =0 truncation, 1: semi coll, 2: Copy closure J+1 = J, P+2 = P)
NL_CLOS = 0;   % nonlinear closure model (0: =0 nmax = jmax, 1: nmax = jmax-j, >1 : nmax = NL_CLOS)
KERN    = 0;   % Kernel model (0 : GK)
INIT_PHI= 0;   % Start simulation with a noisy phi
INIT_ZF = 0;   % Start simulation with a noisy phi
%% OUTPUTS
W_DOUBLE = 0;
W_GAMMA  = 0;
W_PHI    = 0;
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
%% PARAMETER SCANS

if 0
%% Parameter scan over PJ
% PA = [2 4 6 10];
% JA = [1 2 3  5];
PA = [6];
JA = [3];
DTA= DT*ones(size(JA));%./sqrt(JA);
% DTA= DT;
mup_ = MU_P;
muj_ = MU_J;
Nparam = numel(PA);
param_name = 'PJ';
gamma_Ni00 = zeros(Nparam,floor(N/2)+1);
gamma_Nipj = zeros(Nparam,floor(N/2)+1);
Bohm_transport = zeros(Nparam,1);
Ni00_ST  = zeros(Nparam,floor(N/2)+1,floor(SPS2D*TMAX));
for i = 1:Nparam
    % Change scan parameter
    PMAXE = PA(i); PMAXI = PA(i);
    JMAXE = JA(i); JMAXI = JA(i);
    DT = DTA(i);
    setup
    % Run linear simulation
%     system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 1 ./../../../bin/helaz_3 1 1; cd ../../../wk'])
    system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 6 ./../../../bin/helaz_3 1 6; cd ../../../wk'])
%     system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 6 ./../../../bin/helaz_3 2 3; cd ../../../wk'])
%     Load and process results
    %%
    filename = ['../results/',SIMID,'/',PARAMS,'/outputs_00.h5'];
    load_results
    tend   = Ts3D(end); tstart   = 0.4*tend;
    [~,itstart] = min(abs(Ts3D-tstart));
    [~,itend]   = min(abs(Ts3D-tend));
    for ikx = 1:N/2+1
        gamma_Ni00(i,ikx) = (LinearFit_s(Ts3D(itstart:itend)',(squeeze(abs(Ni00(ikx,1,itstart:itend))))));
        Ni00_ST(i,ikx,1:numel(Ts3D)) = squeeze((Ni00(ikx,1,:)));
    end
    tend   = Ts5D(end); tstart   = 0.4*tend;
    [~,itstart] = min(abs(Ts5D-tstart));
    [~,itend]   = min(abs(Ts5D-tend));
    for ikx = 1:N/2+1
        gamma_Nipj(i,ikx) = LinearFit_s(Ts5D(itstart:itend)',squeeze(max(max(abs(Nipj(:,:,ikx,1,itstart:itend)),[],1),[],2)));
    end
    gamma_Ni00(i,:) = real(gamma_Ni00(i,:) .* (gamma_Ni00(i,:)>=0.0));
    gamma_Nipj(i,:) = real(gamma_Nipj(i,:) .* (gamma_Nipj(i,:)>=0.0));
%     kymax = abs(kx(ikymax));
%     Bohm_transport(i) = ETAB/ETAN*gmax/kymax^2;
    % Clean output
    system(['rm -r ',BASIC.RESDIR]);
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
        semilogx(plt(SCALE*kx(2:numel(kx))),plt(gamma_Ni00(i,2:end)),...
            'Color',clr,...
            'LineStyle',linestyle{1},'Marker','^',...
            'DisplayName',['$\eta=',num2str(ETAB/ETAN),'$, $\nu_{',CONAME,'}=',num2str(NU),'$, $P=',num2str(PA(i)),'$, $J=',num2str(JA(i)),'$']);
        hold on;
    end
    grid on; xlabel('$k_z\rho_s^{R}$'); ylabel('$\gamma(N_i^{00})L_\perp/c_s$'); xlim([0.0,max(kx)]);
    title(['$\eta=',num2str(ETAB/ETAN),'$, $\nu_{',CONAME,'}=',num2str(NU),'$'])
    legend('show'); xlim([0.01,10])
saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.fig']);
saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.png']);
end
end

if 1
%% Parameter scan over CO
P=2; J=1;
N       = 20;     % Frequency gridpoints (Nkx = N/2)
L       = 75;     % Size of the squared frequency domain
TMAX  = 200;
DT    = 0.01;
CO_A = [4];
CONAME_A = {};
Nparam = numel(CO_A);
param_name = 'CO';
ETAN    = 2.0;
NU      = 1e-1;   % Collision frequency
PMAXE = P; PMAXI = P;
JMAXE = J; JMAXI = J;
Bohm_transport = zeros(Nparam,1);
gamma_Ni00 = zeros(Nparam,floor(N/2)+1);
gamma_Nipj = zeros(Nparam,floor(N/2)+1);

for i = 1:Nparam
    % Change scan parameter
    CO = CO_A(i);
    setup
    CONAME_A{i} = CONAME;
    % Run linear simulation
    system(...
        ['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 6 ./../../../bin/helaz_3 1 6; cd ../../../wk']...
    )
    %% Load an process results
    filename = ['../results/',SIMID,'/',PARAMS,'/outputs_00.h5'];
    load_results
    tend   = Ts3D(end); tstart   = 0.4*tend;
    [~,itstart] = min(abs(Ts3D-tstart));
    [~,itend]   = min(abs(Ts3D-tend));
    for ikx = 1:N/2+1
        gamma_Ni00(i,ikx) = (LinearFit_s(Ts3D(itstart:itend)',(squeeze(abs(Ni00(ikx,1,itstart:itend))))));
        Ni00_ST(i,ikx,1:numel(Ts3D)) = squeeze((Ni00(ikx,1,:)));
    end
    tend   = Ts5D(end); tstart   = 0.4*tend;
    [~,itstart] = min(abs(Ts5D-tstart));
    [~,itend]   = min(abs(Ts5D-tend));
    for ikx = 1:N/2+1
        gamma_Nipj(i,ikx) = LinearFit_s(Ts5D(itstart:itend)',squeeze(max(max(abs(Nipj(:,:,ikx,1,itstart:itend)),[],1),[],2)));
    end
    gamma_Ni00(i,:) = real(gamma_Ni00(i,:) .* (gamma_Ni00(i,:)>=0.0));
    gamma_Nipj(i,:) = real(gamma_Nipj(i,:) .* (gamma_Nipj(i,:)>=0.0));
%     kymax = abs(kx(ikymax));
%     Bohm_transport(i) = ETAB/ETAN*gmax/kymax^2;
    % Clean output
    system(['rm -r ',BASIC.RESDIR]);
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
        semilogx(plt(SCALE*kx(2:numel(kx))),plt(gamma_Ni00(i,2:end)),...
            'Color',clr,...
            'LineStyle',linestyle{1},'Marker','^',...
            'DisplayName',[CONAME_A{i}]);
        hold on;
    end
    grid on; xlabel('$k_z\rho_s^{R}$'); ylabel('$\gamma(N_i^{00})L_\perp/c_s$'); xlim([0.0,max(kx)]);
    title(['$\eta=',num2str(ETAB/ETAN),'$, $\nu=',num2str(NU),'$',', $P=',num2str(PMAXI),'$, $J=',num2str(JMAXI),'$'])
    legend('show'); xlim([0.01,10])
saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.fig']);
saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.png']);
end

if 0
%% Plot
fig = figure; FIGNAME = 'mixing_length';
plot(eta_B, Bohm_transport)
grid on; xlabel('$L_n/L_B$'); ylabel('$\eta\gamma_{max}/k_{max}^2$');
title(['$P_e=',num2str(PMAXE),'$',', $J_e=',num2str(JMAXE),'$',...
       ', $P_i=',num2str(PMAXE),'$',', $J_i=',num2str(JMAXI),'$'])
saveas(fig,[SIMDIR,FIGNAME,'_vs_',param_name,'_',PARAMS,'.fig']);
end
%%
end