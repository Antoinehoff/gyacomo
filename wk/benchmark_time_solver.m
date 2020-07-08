clear all; close all;
SIMID = 'test_full_coulomb'; % Name of the simulations
addpath(genpath('../matlab')) % ... add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% outputs options
OUTPUTS.nsave_0d = 0;
OUTPUTS.nsave_1d = 0;
OUTPUTS.nsave_2d = 1;
OUTPUTS.nsave_5d = 1;
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
GRID.nkz   = 1;
GRID.kzmin = 1.0;
GRID.kzmax = 1.0;
%% Model parameters
MODEL.CO      = -1;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
MODEL.nu      = 1.0; % collisionality nu*L_perp/Cs0
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
MODEL.eta_T   = 0.0;        % Temperature
MODEL.eta_B   = 0.0;        % Magnetic
% Debye length
MODEL.lambdaD = 0.0;
%% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = '''RK4''';
BASIC.nrun                = 100000;
BASIC.dt                  = 0.05;
BASIC.tmax                = 10.0;
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
run ../matlab/MOLI_time_solver
if params.RK4; MOLIsolvername = 'RK4';
else;          MOLIsolvername = 'ode15i';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis and basic figures
default_plots_options
SAVEFIG = 1;
filename = 'results_00.h5';
default_plots_options
TITLE  = [TITLE,', $k_z=',num2str(GRID.kzmin),'$'];

bare = @(p_,j_) (GRID.jmaxe+1)*p_ + j_ + 1;
bari = @(p_,j_) bare(GRID.pmaxe, GRID.jmaxe) + (GRID.jmaxi+1)*p_ + j_ + 1;
%% Load moments

moment = 'moments_i';

kr     = h5read(filename,['/data/var5d/' moment '/coordkr']);
kz     = h5read(filename,['/data/var5d/' moment '/coordkz']);
time   = h5read(filename,'/data/var5d/time');
Nipj   = zeros(GRID.pmaxi+1, GRID.jmaxi+1,numel(kr),numel(kz),numel(time));
for it = 1:numel(time)
    tmp          = h5read(filename,['/data/var5d/', moment,'/', num2str(it,'%06d')]);
    Nipj(:,:,:,:,it) = tmp.real + 1i * tmp.imaginary;
end

moment = 'moments_e';

kr     = h5read(filename,['/data/var5d/' moment '/coordkr']);
kz     = h5read(filename,['/data/var5d/' moment '/coordkz']);
time   = h5read(filename,'/data/var5d/time');
Nepj   = zeros(GRID.pmaxe+1, GRID.jmaxe+1,numel(kr),numel(kz),numel(time));
for it = 1:numel(time)
    tmp          = h5read(filename,['/data/var5d/', moment,'/', num2str(it,'%06d')]);
    Nepj(:,:,:,:,it) = tmp.real + 1i * tmp.imaginary;
end


%% Plot moments
fig = figure;

x1 = time;
x2 = results.time;
ic = 1;

% Electrons
subplot(321)
for ip = 0:1
    for ij = 0:1
        y1 = squeeze(real(Nepj(ip+1,ij+1,1,1,:)));
        plot(x1,y1,'-','DisplayName',...
            ['HeLaZ $N_e^{',num2str(ip),num2str(ij),'}$'],...
            'color', line_colors(ic,:))
        hold on
        y2 = squeeze(real(results.Napj(:,bare(ip,ij))));
        plot(x2,y2,'--','DisplayName',...
            ['MOLI $N_e^{',num2str(ip),num2str(ij),'}$'],...
            'color', line_colors(ic,:))
        hold on
        ic = ic + 1;
    end
end
grid on
xlabel('$t$')
ylabel(['Re$(N_e^{pj})$'])

% Ions
ic = 1;
subplot(322)
for ip = 0:1
    for ij = 0:1
        y1 = squeeze(real(Nipj(ip+1,ij+1,1,1,:)));
        plot(x1,y1,'-','DisplayName',...
            ['HeLaZ $N_i^{',num2str(ip),num2str(ij),'}$'],...
            'color', line_colors(ic,:))
        hold on
        y2 = squeeze(real(results.Napj(:,bari(ip,ij))));
        plot(x2,y2,'--','DisplayName',...
            ['MOLI $N_i^{',num2str(ip),num2str(ij),'}$'],...
            'color', line_colors(ic,:))
        hold on
        ic = ic + 1;
    end
end
grid on
xlabel('$t$')
ylabel(['Re$(N_i^{pj})$'])
%suptitle(TITLE);
%if SAVEFIG; FIGNAME = ['Nipj_kz_',num2str(GRID.kzmin)]; save_figure; end;

% phi
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

%fig = figure;
%Real
subplot(323)
plot(timephi,real(phiHeLaZ),'-','DisplayName','HeLaZ RK4')
hold on
plot(timephiMOLI,real(phiMOLI),'--','DisplayName',['MOLI ',MOLIsolvername])
grid on
xlabel('$t$')
ylabel('Re$(\phi)$')
%Imag
subplot(324)
plot(timephi,imag(phiHeLaZ),'-','DisplayName','HeLaZ RK4')
hold on
plot(timephiMOLI,imag(phiMOLI),'--','DisplayName',['MOLI ',MOLIsolvername])
grid on
xlabel('$t$')
ylabel('Im$(\phi)$')
%if SAVEFIG; FIGNAME = ['phi_kz_',num2str(GRID.kzmin)]; save_figure; end;

%% phi error
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

ERR1 = abs(real(phiMOLI(1:numel(timephi)) - phiHeLaZ));
ERR2 = abs(imag(phiMOLI(1:numel(timephi)) - phiHeLaZ));

%fig = figure;
subplot(325);
plot(timephi,ERR1,'-','DisplayName','Real')
hold on
plot(timephi,ERR2,'--','DisplayName','Imag')
%title(TITLE);
grid on
xlabel('$t$')
ylabel('$e(\phi)$')
legend('show')
%if SAVEFIG; FIGNAME = ['err_phi_kz_',num2str(GRID.kzmin)]; save_figure; end;

% Growth rate fit quantity (Ni00)
subplot(326);

y1 = squeeze(abs(Nipj(1,1,1,1,:))); % HeLaZ
y2 = squeeze(abs(results.Napj(:,bari(0,0))));
ERR3  = abs(y2(1:numel(timephi))-y1);

yyaxis left
semilogy(x1,y1,'-','DisplayName','HeLaZ',...
    'color', line_colors(1,:)); hold on;
semilogy(x2,y2,'--','DisplayName','MOLI',...
    'color', 'k')
grid on
xlabel('$t$')
ylabel('$|N_i^{00}|$')
legend('show')

yyaxis right
semilogy(x1, ERR3,'color', line_colors(2,:));
grid on
xlabel('$t$')
ylabel('$e(N_i^{00})$')

% Finish and save
suptitle(TITLE);
if SAVEFIG; FIGNAME = ['kz_',num2str(GRID.kzmin)]; save_figure; end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
