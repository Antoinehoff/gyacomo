clear all;
addpath(genpath('../matlab')) % ... add
default_plots_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 1e-2;   % Collision frequency
TAU     = 1.0;    % e/i temperature ratio
ETAB    = 1.0;
ETAN    = 1.0;    % Density gradient
ETAT    = 0.0;    % Temperature gradient
MU      = 1e-4;   % Hyper diffusivity coefficient
LAMBDAD = 0.0; 
NOISE0  = 1.0e-5;
%% GRID PARAMETERS
N       = 50;     % Frequency gridpoints (Nkr = N/2)
L       = 30;     % Size of the squared frequency domain
PMAXE   = 20;
JMAXE   = 10; 
PMAXI   = 20;
JMAXI   = 10;
KREQ0   = 1;      % put kr = 0
KPAR    = 0.0;    % Parellel wave vector component
%% TIME PARAMETERS 
TMAX    = 150;  % Maximal time unit
DT      = 5e-3;   % Time step
SPS2D   = 1;      % Sampling per time unit for 2D arrays
SPS5D   = 0;    % Sampling per time unit for 5D arrays
RESTART = 0;      % To restart from last checkpoint
JOB2LOAD= 00;
%% OPTIONS
SIMID   = 'linear_study';  % Name of the simulation
NON_LIN = 1 *(1-KREQ0);   % activate non-linearity (is cancelled if KREQ0 = 1)
CO      = -2;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unused
KR0KH   = 0; A0KH = 0; % Background phi mode to drive Ray-Tay inst.
NO_E    = 0;  % Remove electrons dynamic
% DK    = 0;  % Drift kinetic model (put every kernel_n to 0 except n=0 to 1)
JOBNUM = 00;

setup


%% PARAMETER SCANS
if 0
%% Parameter scan over PJ
ETAB  = 1.0;
TMAX  = 150;
PA = [2,3,4,5,6,8,10,12,14,16,20];
JA = [1,1,2,2,3,4,5,6,7,8,10];
% PA = [40];
% JA = [20];
Nparam = numel(PA);
param_name = 'PJ';

gamma_Ni = zeros(Nparam,N);

for i = 1:Nparam
    % Change scan parameter
    PMAXE = PA(i); PMAXI = PA(i);
    JMAXE = JA(i); JMAXI = JA(i);
    setup
    % Run linear simulation
    run
    % Load and process results
    load_results
    tend   = Ts2D(end); tstart   = 0.6*tend; 
    for ikz = 1:N
        gamma_Ni(i,ikz) = LinearFit_s(Ts2D,squeeze(abs(Ni00(1,ikz,:))),tstart,tend);
    end
    gamma_Ni(i,:) = real(gamma_Ni(i,:) .* (gamma_Ni(i,:)>=0.0));
    % Clean output
    system(['rm -r ',BASIC.RESDIR])
end

if 1
%% Plot
fig = figure; FIGNAME = 'linear_study';
plt = @(x) circshift(x,N/2-1);
for i = 1:Nparam
    clr       = line_colors(mod(i-1,numel(line_colors(:,1)))+1,:);
    linestyle = line_styles(floor((i-1)/numel(line_colors(:,1)))+1);
    plot(plt(kz),plt(gamma_Ni(i,:)),...
        'Color',clr,...
        'LineStyle',linestyle{1},...
        'DisplayName',['$P=$',num2str(PA(i)),', $J=$',num2str(JA(i))]);
    hold on;
end
grid on; xlabel('$k_z$'); ylabel('$\gamma(N_i^{00})$'); xlim([0.0,max(kz)]);
title(['$\eta_B=',num2str(ETAB),'$'])
legend('show')
saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.fig']);
end
end
if 0
%% Parameter scan over eta_B
PMAXE = 20; PMAXI = 20;
JMAXE = 10; JMAXI = 10;
TMAX  = 50;
eta_B = [0.1, 0.33, 0.5, 0.67, 1.0];
Nparam = numel(eta_B);
param_name = 'etaB';

gamma_Ni = zeros(Nparam,N);

for i = 1:Nparam
    % Change scan parameter
    ETAB = eta_B(i);
    setup
    % Run linear simulation
    run
    % Load and process results
    load_results
    tend   = Ts2D(end); tstart   = 0.6*tend; 
    for ikz = 1:N
        gamma_Ni(i,ikz) = LinearFit_s(Ts2D,squeeze(abs(Ni00(1,ikz,:))),tstart,tend);
    end
    gamma_Ni(i,:) = real(gamma_Ni(i,:) .* (gamma_Ni(i,:)>=0.0));
    % Clean output
    system(['rm -r ',BASIC.RESDIR])
end

if 1
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
grid on; xlabel('$k_z$'); ylabel('$\gamma(N_i^{00})$'); xlim([0.0,max(kz)]);
title(['$P_e=',num2str(PMAXE),'$',', $J_e=',num2str(JMAXE),'$',...
       ', $P_i=',num2str(PMAXE),'$',', $J_i=',num2str(JMAXI),'$'])
legend('show')
saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.fig']);
end
end