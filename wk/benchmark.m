SIMID = 'benchmark'; % Name of the simulations
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
GRID.pmaxe = 20;
GRID.jmaxe = 10;
GRID.pmaxi = 20;
GRID.jmaxi = 10;
GRID.nkr   = 1;
GRID.krmin = 0.;
GRID.krmax = 0.;
GRID.nkz   = 1;
GRID.kzmin = 1.0;
GRID.kzmax = 1.0;
%% Model parameters
MODEL.CO      = -1;  % Collision operator (-1 = Full Coulomb, 0 = Dougherty)
MODEL.nu      = 0.0; % collisionality nu*L_perp/Cs0
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
MODEL.eta_n   = 0.0;        % Density
MODEL.eta_T   = 1.0;        % Temperature
MODEL.eta_B   = 0.0;        % Magnetic
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
run ../matlab/MOLI_time_solver
if params.RK4; MOLIsolvername = 'RK4';
else;          MOLIsolvername = 'ode15i';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis and basic figures
SAVEFIG = 1;
filename = 'results_00.h5';
TITLE  = [];
TITLE = [TITLE,'$\eta_n=',num2str(MODEL.eta_n),'$, '];
TITLE = [TITLE,'$\eta_B=',num2str(MODEL.eta_B),'$, '];
TITLE = [TITLE,'$\eta_T=',num2str(MODEL.eta_T),'$, '];
TITLE = [TITLE,   '$\nu=',num2str(MODEL.nu),'$, '];
TITLE = [TITLE, '$(P,J)=(',num2str(GRID.pmaxe),',',num2str(GRID.jmaxe),')$'];
%% Nipj

moment = 'Ni00';

kr       = h5read(filename,['/data/var2d/' moment '/coordkr']);
kz       = h5read(filename,['/data/var2d/' moment '/coordkz']);
timeNi     = h5read(filename,'/data/var2d/time');
Nipj     = zeros(numel(timeNi),numel(kr),numel(kz));
for it = 1:numel(timeNi)
    tmp          = h5read(filename,['/data/var2d/', moment,'/', num2str(it,'%06d')]);
    Nipj(it,:,:) = tmp.real + 1i * tmp.imaginary; 
end

%% phi
timephi  = h5read(filename,'/data/var2d/time');
kr       = h5read(filename,'/data/var2d/phi/coordkr');
kz       = h5read(filename,'/data/var2d/phi/coordkz');
phiHeLaZ      = zeros(numel(timephi),numel(kr),numel(kz));
for it = 1:numel(timephi)
    tmp         = h5read(filename,['/data/var2d/phi/' num2str(it,'%06d')]);
    phiHeLaZ(it,:,:) = tmp.real + 1i * tmp.imaginary;
end

timephiMOLI = results.time;
phiMOLI  = zeros(size(timephiMOLI));
for it = 1:numel(timephiMOLI)
    phiMOLI(it) = get_phi(results.Napj(it,:),params,options);
end

fig = figure;
%HeLaZ results
semilogy(timephi,abs(phiHeLaZ),'-','DisplayName','HeLaZ RK4')
hold on
title(TITLE);
%MOLI results
semilogy(timephiMOLI,abs(phiMOLI),'--','DisplayName',['MOLI ',MOLIsolvername])
grid on
xlabel('$t$')
ylabel('$|\phi|$')
legend('show')
if SAVEFIG; FIGNAME = 'phi'; save_figure; end;


fig = figure;
%HeLaZ results
x1 = timeNi;
y1 = abs(Nipj);
semilogy(x1,y1,'-','DisplayName','HeLaZ RK4')
hold on
%MOLI results
x2 = results.time;
y2 = abs(results.Napj(:,1));
semilogy(x2,y2,'--','DisplayName',['MOLI ',MOLIsolvername]);
title(TITLE);
grid on
legend('show')
xlabel('$t$')
ylabel(['$|',moment,'|$'])
if SAVEFIG; FIGNAME = 'Ni00'; save_figure; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%