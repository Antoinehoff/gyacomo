%clear all;
addpath(genpath('../matlab')) % ... add
default_plots_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 1e-1;   % Collision frequency
TAU     = 1.0;    % e/i temperature ratio
ETAB    = 0.66;
ETAN    = 1.0;    % Density gradient
ETAT    = 0.0;    % Temperature gradient
MU      = 1e-3;   % Hyper diffusivity coefficient
LAMBDAD = 0.0;
NOISE0  = 1.0e-5;
%% GRID PARAMETERS
N       = 100;     % Frequency gridpoints (Nkr = N/2)
L       = 100;     % Size of the squared frequency domain
KREQ0   = 1;      % put kr = 0
%% TIME PARAMETERS
TMAX    = 100;  % Maximal time unit
DT      = 1e-2;   % Time step
SPS0D   = 0.5;      % Sampling per time unit for 2D arrays
SPS2D   = 1;      % Sampling per time unit for 2D arrays
SPS5D   = 0.1;    % Sampling per time unit for 5D arrays
SPSCP   = 0;    % Sampling per time unit for checkpoints
RESTART = 0;      % To restart from last checkpoint
JOB2LOAD= 00;
%% OPTIONS
SIMID   = 'linear_study';  % Name of the simulation
NON_LIN = 0 *(1-KREQ0);   % activate non-linearity (is cancelled if KREQ0 = 1)
CO      = -1;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unused
NO_E    = 0;  % Remove electrons dynamic
% DK    = 0;  % Drift kinetic model (put every kernel_n to 0 except n=0 to 1)
JOBNUM = 00;
KPAR    = 0.0;    % Parellel wave vector component

%% PARAMETER SCANS
if 1
%% Parameter scan over PJ
PA = [2, 6, 8, 12];
JA = [1, 3, 4,  6];
% PA = [2];
% JA = [1];
Nparam = numel(PA);

gamma_Ni = zeros(Nparam,N/2+1);
Ni00_ST  = zeros(Nparam,N/2+1,TMAX);
for i = 1:Nparam
    % Change scan parameter
    PMAXE = PA(i); PMAXI = PA(i);
    JMAXE = JA(i); JMAXI = JA(i);
    setup
    % Run linear simulation
    system(...
        ['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 8 ./../../../bin/helaz; cd ../../../wk']...
    )
    % Load and process results
    load_results
    tend   = Ts2D(end); tstart   = 0.4*tend;
    for ikr = 1:N/2+1
        gamma_Ni(i,ikr) = LinearFit_s(Ts2D,squeeze(abs(Ni00(ikr,1,:))),tstart,tend);
        Ni00_ST(i,ikr,1:numel(Ts2D)) = squeeze((Ni00(ikr,1,:)));
    end
    gamma_Ni(i,:) = real(gamma_Ni(i,:) .* (gamma_Ni(i,:)>=0.0));
    % Clean output
    %system(['rm -r ',BASIC.RESDIR])
end
% if 0
% %% Plot
% fig = figure; FIGNAME = 'space_time_Ni00';
% i = 1;
%
% [YY, XX] = meshgrid(Ts2D,kz);
% plt = @(x) squeeze(log(abs(x)))
% pclr = pcolor(XX,YY,plt(Ni00_ST(i,:,:)));set(pclr, 'edgecolor','none');
% grid on; xlabel('$k_z$'); ylabel('$t$');
% title(['$\eta_B=',num2str(ETAB),'$,$P,J=',num2str(PA(i)),',',JA(i),'$'])
% legend('show')
% saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.fig']);
% saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.png']);
% end
if 0
%% Plot
fig = figure; FIGNAME = 'space_time_Ni00';
i = 3;
plt = @(x) squeeze((abs(x)));
for ikr=1:4:numel(kr)/2+1
    NAME = ['$k_z=',num2str(kr(ikr)),'$'];
    semilogy(1:TMAX,plt(Ni00_ST(i,ikr,1:TMAX)), 'DisplayName', NAME); hold on;
end
grid on; xlabel('$t c_s/\rho_s$'); ylabel('$|Ni^{00}|$');
title(['$\eta_B=',num2str(ETAB),'$; $P,J=',num2str(PA(i)),',',num2str(JA(i)),'$, $\nu=',num2str(NU),'$'])
legend('show');
saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.fig']);
saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.png']);
end
if 1
%% Plot
fig = figure; FIGNAME = 'linear_study';
plt = @(x) x;
for i = 1:Nparam
    clr       = line_colors(mod(i-1,numel(line_colors(:,1)))+1,:);
    linestyle = line_styles(floor((i-1)/numel(line_colors(:,1)))+1);
    plot(plt(kr),plt(gamma_Ni(i,:)),...
        'Color',clr,...
        'LineStyle',linestyle{1},...
        'DisplayName',['$P=$',num2str(PA(i)),', $J=$',num2str(JA(i))]);
    hold on;
end
grid on; xlabel('$k_z$'); ylabel('$\gamma(N_i^{00})$'); xlim([0.0,max(kr)]);
title(['$\eta_B=',num2str(ETAB),'$, $\nu=',num2str(NU),'$'])
legend('show')
saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.fig']);
saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.png']);
end
end
if 0
%% Parameter scan over eta_B
PMAXE = 20; PMAXI = 20;
JMAXE = 10; JMAXI = 10;
TMAX  = 200;
eta_B = [0.1, 0.33, 0.5, 0.67, 1.0];
Nparam = numel(eta_B);
param_name = 'etaB';
CO      = -2;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
NU      = 1e-2;   % Collision frequency

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
grid on; xlabel('$k_z\rho_s$'); ylabel('$\gamma(N_i^{00})\rho_2/c_s$'); xlim([0.0,max(kz)]);
title(['$P_e=',num2str(PMAXE),'$',', $J_e=',num2str(JMAXE),'$',...
       ', $P_i=',num2str(PMAXE),'$',', $J_i=',num2str(JMAXI),'$'])
legend('show')
saveas(fig,[SIMDIR,'gamma_Ni_vs_',param_name,'_',PARAMS,'.fig']);
end
end
