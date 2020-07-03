SIMID = 'benchmark_kperp_scan'; % Name of the simulations
addpath(genpath('../matlab')) % ... add 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% outputs options
OUTPUTS.nsave_0d = 0;
OUTPUTS.nsave_1d = 0;
OUTPUTS.nsave_2d = 1;
OUTPUTS.nsave_5d = 0;
OUTPUTS.write_Ni00    = '.true.';
OUTPUTS.write_moments = '.true.';
OUTPUTS.write_phi     = '.true.';
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
MODEL.CO      = -2;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
MODEL.nu      = 0.1; % collisionality nu*L_perp/Cs0
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
BASIC.tmax                = 100.0;
INITIAL.initback_moments  = 0.01;
INITIAL.initnoise_moments = 0.;
INITIAL.iseed             = 42;

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

if     MODEL.CO == -1; CONAME = 'FC';
elseif MODEL.CO == -2; CONAME = 'DC';
elseif MODEL.CO ==  0; CONAME = 'LB'; end;

TITLE  = [];
TITLE = [TITLE,'$\eta_n=',num2str(MODEL.eta_n),'$, '];
TITLE = [TITLE,'$\eta_B=',num2str(MODEL.eta_B),'$, '];
TITLE = [TITLE,'$\eta_T=',num2str(MODEL.eta_T),'$, '];
TITLE = [TITLE,   '$\nu=',num2str(MODEL.nu),'$, '];
TITLE = [TITLE, '$(P,J)=(',num2str(GRID.pmaxe),',',num2str(GRID.jmaxe),')$'];

%% Growth rate analysis
gammas = zeros(numel(kr),numel(kz));
shifts = zeros(numel(kr),numel(kz));

moment = 'Ni00';

kr       = h5read(filename,['/data/var2d/' moment '/coordkr']);
kz       = h5read(filename,['/data/var2d/' moment '/coordkz']);
timeNi   = h5read(filename,'/data/var2d/time');
Nipj     = zeros(numel(timeNi),numel(kr),numel(kz));

for it = 1:numel(timeNi)
    tmp          = h5read(filename,['/data/var2d/', moment,'/', num2str(it,'%06d')]);
    Nipj(it,:,:) = tmp.real + 1i * tmp.imaginary; 
end

% Linear fit of log(Napj)
x1    = timeNi;
itmin = ceil(0.5 * numel(timeNi)); %Take the second half of the time evolution
ikr   = 1;

for ikz = 1:numel(kz)
    fit = polyfit(x1(itmin:end),log(abs(Nipj(itmin:end,ikr,ikz))),1);
    gammas(ikr,ikz) = fit(1);
    shifts(ikr,ikz) = fit(2);
end

%% Plot
fig = figure;

%HeLaZ results
X = kz*sqrt(1+MODEL.tau_i);
Y = gammas(1,:);
plot(X,Y,'-','DisplayName','HeLaZ')
hold on

%MOLI results
X = kz*sqrt(1+MODEL.tau_i);
Y = results.kperp.Maxgammas;
plot(X,Y,'--','DisplayName','MOLI');
title(TITLE);
grid on
legend('show')
xlabel('$k_\perp * (1+\tau)^{1/2}$')
ylabel('$\gamma L_\perp/c_{s} $')
if SAVEFIG; FIGNAME = 'g_vs_kz'; save_figure; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%