%clear all;
addpath(genpath('../matlab')) % ... add
default_plots_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 1.0;   % Collision frequency
TAU     = 1.0;    % e/i temperature ratio
ETAB    = 0.6;
ETAN    = 1.0;    % Density gradient
ETAT    = 0.0;    % Temperature gradient
NU_HYP  = 0.1;   % Hyperdiffusivity coefficient
LAMBDAD = 0.0;
NOISE0  = 1.0e-5;
%% GRID PARAMETERS
N       = 150;     % Frequency gridpoints (Nkr = N/2)
L       = 70;     % Size of the squared frequency domain
KREQ0   = 1;      % put kr = 0
MU_P    = 0.0;     % Hermite  hyperdiffusivity -mu_p*(d/dvpar)^4 f
MU_J    = 0.0;     % Laguerre hyperdiffusivity -mu_j*(d/dvperp)^4 f
%% TIME PARMETERS
TMAX    = 100;  % Maximal time unit
DT      = 1e-2;   % Time step
SPS0D   = 0.5;      % Sampling per time unit for 2D arrays
SPS2D   = 1;      % Sampling per time unit for 2D arrays
SPS5D   = 1;    % Sampling per time unit for 5D arrays
SPSCP   = 0;    % Sampling per time unit for checkpoints
RESTART = 0;      % To restart from last checkpoint
JOB2LOAD= 00;
%% OPTIONS
SIMID   = 'linear_study_test_mu_kin';  % Name of the simulation
NON_LIN = 0 *(1-KREQ0);   % activate non-linearity (is cancelled if KREQ0 = 1)
CO      = -3;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
CLOS    = 0;   % Closure model (0: =0 truncation, 1: semi coll, 2: Copy closure J+1 = J, P+2 = P)
KERN    = 0;   % Kernel model (0 : GK)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unused
% DK    = 0;  % Drift kinetic model (put every kernel_n to 0 except n=0 to 1)
JOBNUM = 00;
KPAR    = 0.0;    % Parellel wave vector component
HD_CO   = 0.5;    % Hyper diffusivity cutoff ratio
kmax    = N*pi/L;% Highest fourier mode
MU      = NU_HYP/(HD_CO*kmax)^4 % Hyperdiffusivity coefficient

%% PARAMETER SCANS
if 1
%% Parameter scan over PJ
PA = [2, 3, 4, 6, 8, 10];
JA = [1, 2, 2, 3, 4,  5];
DTA= DT./sqrt(JA);
mup_ = MU_P;
muj_ = MU_J;
% PA = [4];
% JA = [2];
Nparam = numel(PA);
param_name = 'PJ';
gamma_Ni00 = zeros(Nparam,N/2+1);
gamma_Ni21 = zeros(Nparam,N/2+1);
Bohm_transport = zeros(Nparam,1);
Ni00_ST  = zeros(Nparam,N/2+1,TMAX);
for i = 1:Nparam
    % Change scan parameter
    PMAXE = PA(i); PMAXI = PA(i);
    JMAXE = JA(i); JMAXI = JA(i);
    DT = DTA(i);
    MU_P = mup_/PMAXE^2;
    MU_J = muj_/JMAXE^3;
    setup
    % Run linear simulation
    system(...
        ['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 4 ./../../../bin/helaz; cd ../../../wk']...
    )
    % Load and process results
    load_results
    tend   = Ts2D(end); tstart   = 0.4*tend;
    for ikr = 1:N/2+1
        gamma_Ni00(i,ikr) = LinearFit_s(Ts2D,squeeze(abs(Ni00(ikr,1,:))),tstart,tend);
        Ni00_ST(i,ikr,1:numel(Ts2D)) = squeeze((Ni00(ikr,1,:)));
    end
    tend   = Ts5D(end); tstart   = 0.4*tend;
    for ikr = 1:N/2+1
        gamma_Ni21(i,ikr) = LinearFit_s(Ts5D,squeeze(abs(Nipj(3,2,ikr,1,:))),tstart,tend);
    end
    gamma_Ni00(i,:) = real(gamma_Ni00(i,:) .* (gamma_Ni00(i,:)>=0.0));
    gamma_Ni21(i,:) = real(gamma_Ni21(i,:) .* (gamma_Ni21(i,:)>=0.0));
    kzmax = abs(kr(ikzmax));
    Bohm_transport(i) = ETAB/ETAN*gmax/kzmax^2;
    % Clean output
    system(['rm -r ',BASIC.RESDIR])
end

if 1
%% Plot
fig = figure; FIGNAME = 'linear_study';
plt = @(x) x;
subplot(211)
    for i = 1:Nparam
        clr       = line_colors(mod(i-1,numel(line_colors(:,1)))+1,:);
        linestyle = line_styles(floor((i-1)/numel(line_colors(:,1)))+1);
        plot(plt(kr),plt(gamma_Ni00(i,:)),...
            'Color',clr,...
            'LineStyle',linestyle{1},...
            'DisplayName',['$P=$',num2str(PA(i)),', $J=$',num2str(JA(i))]);
        hold on;
    end
    grid on; xlabel('$k_z\rho_s$'); ylabel('$\gamma(N_i^{00})\rho_s/c_s$'); xlim([0.0,max(kr)]);
    title(['$\eta_B=',num2str(ETAB),'$, $\nu_{',CONAME,'}=',num2str(NU),'$, ', CLOSNAME])
    legend('show')
subplot(212)
    for i = 1:Nparam
        clr       = line_colors(mod(i-1,numel(line_colors(:,1)))+1,:);
        linestyle = line_styles(floor((i-1)/numel(line_colors(:,1)))+1);
        plot(plt(kr),plt(gamma_Ni21(i,:)),...
            'Color',clr,...
            'LineStyle',linestyle{1},...
            'DisplayName',['$P=$',num2str(PA(i)),', $J=$',num2str(JA(i))]);
        hold on;
    end
    grid on; xlabel('$k_z\rho_s$'); ylabel('$\gamma(N_i^{21})\rho_s/c_s$'); xlim([0.0,max(kr)]);
    title(['$\eta_B=',num2str(ETAB),'$, $\nu_{',CONAME,'}=',num2str(NU),'$, ', CLOSNAME])
    legend('show')
saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.fig']);
saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.png']);
end
end

if 0
%% Parameter scan over eta_B
PMAXE = 6; PMAXI = 6;
JMAXE = 3; JMAXI = 3;
TMAX  = 400;
DT    = 0.001;
eta_B = [0.45, 0.5, 0.55, 0.6, 0.65];
Nparam = numel(eta_B);
param_name = 'etaB';
CO      = -2;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
NU      = 1e-2;   % Collision frequency

Bohm_transport = zeros(Nparam,1);
gamma_Ni = zeros(Nparam,N);

for i = 1:Nparam
    % Change scan parameter
    ETAB = eta_B(i);
    setup
    % Run linear simulation
    system(...
        ['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 4 ./../../../bin/helaz; cd ../../../wk']...
    )
    % Load and process results
    load_results
    tend   = Ts2D(end); tstart   = 0.4*tend;
    for ikr = 1:N/2+1
        gamma_Ni(i,ikr) = LinearFit_s(Ts2D,squeeze(abs(Ni00(ikr,1,:))),tstart,tend);
        Ni00_ST(i,ikr,1:numel(Ts2D)) = squeeze((Ni00(ikr,1,:)));
    end
    gamma_Ni(i,:) = real(gamma_Ni(i,:) .* (gamma_Ni(i,:)>=0.0));
    [gmax,ikzmax] = max(gamma_Ni(i,:));
    kzmax = abs(kr(ikzmax));
    Bohm_transport(i) = ETAB/ETAN*gmax/kzmax^2;
    % Clean output
    system(['rm -r ',BASIC.RESDIR])
end

if 0
%% Plot
fig = figure; FIGNAME = 'linear_study';
plt = @(x) circshift(x,N/2-1);
for i = 1:Nparam
    clr       = line_colors(mod(i-1,numel(line_colors(:,1)))+1,:);
    linestyle = line_styles(floor((i-1)/numel(line_colors(:,1)))+1);
    plot(plt(kz),plt(gamma_Ni(i,:)),...
        'Color',clr,...
        'LineStyle',linestyle{1},...
        'DisplayName',['$\eta_B=$',num2str(eta_B(i))]);
    hold on;
end
grid on; xlabel('$k_z\rho_s$'); ylabel('$\gamma(N_i^{00})\rho_2/c_s$'); xlim([0.0,max(kz)]);
title(['$P_e=',num2str(PMAXE),'$',', $J_e=',num2str(JMAXE),'$',...
       ', $P_i=',num2str(PMAXE),'$',', $J_i=',num2str(JMAXI),'$'])
legend('show')
saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.fig']);
end

if 1
%% Plot
fig = figure; FIGNAME = 'mixing_length';
plot(eta_B, Bohm_transport)
grid on; xlabel('$L_n/L_B$'); ylabel('$\eta\gamma_{max}/k_{max}^2$');
title(['$P_e=',num2str(PMAXE),'$',', $J_e=',num2str(JMAXE),'$',...
       ', $P_i=',num2str(PMAXE),'$',', $J_i=',num2str(JMAXI),'$'])
saveas(fig,[SIMDIR,FIGNAME,'_vs_',param_name,'_',PARAMS,'.fig']);
end
end
