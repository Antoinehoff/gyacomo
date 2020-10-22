SIMID = 'benchmark_kperp_scan'; % Name of the simulations
addpath(genpath('../matlab')) % ... add 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% outputs options
OUTPUTS.nsave_0d = 0;
OUTPUTS.nsave_1d = 0;
OUTPUTS.nsave_2d = 1;
OUTPUTS.nsave_5d = 0;
OUTPUTS.write_Na00    = '.true.';
OUTPUTS.write_moments = '.true.';
OUTPUTS.write_phi     = '.true.';
OUTPUTS.write_non_lin = '.false.';
OUTPUTS.write_doubleprecision = '.true.';
OUTPUTS.resfile0      = '''results''';
%% Grid parameters
GRID.pmaxe = 15;
GRID.jmaxe = 6;
GRID.pmaxi = 15;
GRID.jmaxi = 6;
GRID.nkr   = 1;
GRID.krmin = 0.;
GRID.krmax = 0.;
GRID.nkz   = 20;
GRID.kzmin = 0.1;
GRID.kzmax = 1.5;
%% Model parameters
MODEL.CO      = -1;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
MODEL.NON_LIN = 0;   % Non linear term
MODEL.nu      = 0.01; % collisionality nu*L_perp/Cs0
% temperature ratio T_a/T_e
MODEL.tau_e   = 1.0;
MODEL.tau_i   = 1.0;
% mass ratio sqrt(m_a/m_i)
MODEL.sigma_e = 0.0233380;
MODEL.sigma_i = 1.0;
% charge q_a/e
MODEL.q_e     =-1.0;
MODEL.q_i     = 1.0;
% gradients L_perp/L_x
MODEL.eta_n   = 1.0;        % Density
MODEL.eta_T   = 0.0;        % Temperature
MODEL.eta_B   = 0.5;        % Magnetic
% Debye length
MODEL.lambdaD = 0.0;
%% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = '''RK4''';
BASIC.nrun                = 100000;
BASIC.dt                  = 0.05;
BASIC.tmax                = 50.0;
INITIAL.initback_moments  = 0.01;
INITIAL.initnoise_moments = 0.;
INITIAL.iseed             = 42;
INITIAL.selfmat_file = ...
    ['''../iCa/self_Coll_GKE_0_GKI_0_ESELF_1_ISELF_1_Pmaxe_',num2str(GRID.pmaxe),...
    '_Jmaxe_',num2str(GRID.jmaxe),'_Pmaxi_',num2str(GRID.pmaxi),'_Jmaxi_',...
    num2str(GRID.jmaxi),'_pamaxx_10.h5'''];
INITIAL.eimat_file = ...
    ['''../iCa/ei_Coll_GKE_0_GKI_0_ETEST_1_EBACK_1_Pmaxe_',num2str(GRID.pmaxe),...
    '_Jmaxe_',num2str(GRID.jmaxe),'_Pmaxi_',num2str(GRID.pmaxi),'_Jmaxi_',...
    num2str(GRID.jmaxi),'_pamaxx_10_tau_1.0000_mu_0.0233.h5'''];
INITIAL.iemat_file = ...
    ['''../iCa/ie_Coll_GKE_0_GKI_0_ITEST_1_IBACK_1_Pmaxe_',num2str(GRID.pmaxe),...
    '_Jmaxe_',num2str(GRID.jmaxe),'_Pmaxi_',num2str(GRID.pmaxi),'_Jmaxi_',...
    num2str(GRID.jmaxi),'_pamaxx_10_tau_1.0000_mu_0.0233.h5'''];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write input file
INPUT = write_fort90(OUTPUTS,GRID,MODEL,INITIAL,TIME_INTEGRATION,BASIC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run HeLaZ
nproc = 1;
EXEC  = ' ../bin/helaz ';
RUN   = ['mpirun -np ' num2str(nproc)];
CMD   = [RUN, EXEC, INPUT];
system(CMD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run MOLI with same input parameters
params.RK4 = 1;
run ../matlab/MOLI_kperp_scan
if params.RK4; MOLIsolvername = 'RK4';
else;          MOLIsolvername = 'ode15i';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis and basic figures
SAVEFIG = 1;
filename = 'results_00.h5';
default_plots_options

%% Growth rate analysis
moment = 'Ni00';

kr       = h5read(filename,['/data/var2d/' moment '/coordkr']);
kz       = h5read(filename,['/data/var2d/' moment '/coordkz']);
timeNi   = squeeze(h5read(filename,'/data/var2d/time'));
Ni00     = zeros(numel(kr),numel(kz),numel(timeNi));

for it = 1:numel(timeNi)
    tmp          = h5read(filename,['/data/var2d/', moment,'/', num2str(it,'%06d')]);
    Ni00(:,:,it) = tmp.real + 1i * tmp.imaginary; 
end

% Amplitude ratio method
ikr   = 1;
gammas = zeros(numel(kr),numel(kz));
Napj  = ones(numel(timeNi),(GRID.pmaxe+1)*(GRID.jmaxe+1)+(GRID.pmaxi+1)*(GRID.jmaxi+1));
for ikz = 1:numel(kz)
   gammas(ikr,ikz) = LinearFit_s(squeeze(timeNi), squeeze(abs(Ni00(ikr,ikz,:))));
end

%% Plot
fig = figure;

X = kz*sqrt(1+MODEL.tau_i); % Convert to Ricci 2006 Normalization

subplot(121) % growth rate
%HeLaZ results
Y1 = gammas(1,:);
plot(X,Y1,'-','DisplayName','HeLaZ')
hold on

%MOLI results
Y2 = results.kperp.Maxgammas;
plot(X,Y2,'--','DisplayName','MOLI');
title(TITLE);
grid on
legend('show')
xlabel('$k_\perp * (1+\tau)^{1/2}$')
ylabel('$\gamma L_\perp/c_{s} $')

subplot(122) % Error
ERR = abs(gammas(1,:)' - results.kperp.Maxgammas);

plot(X,ERR,'-','DisplayName','MOLI relative error');
grid on
legend('show')
xlabel('$k_\perp * (1+\tau)^{1/2}$')
ylabel('$\epsilon_\gamma $')

if SAVEFIG; FIGNAME = 'g_vs_kz'; save_figure; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%