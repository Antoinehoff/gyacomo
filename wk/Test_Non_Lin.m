%% Run linear simulation on a kr,kz grid and compute the non linear term afterwards
clear all; close all;
BASIC.SIMID = 'test_non_lin'; % Name of the simulations
addpath(genpath('../matlab')) % ... add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% outputs options
OUTPUTS.nsave_0d = 0;
OUTPUTS.nsave_1d = 0;
OUTPUTS.nsave_2d = 1;
OUTPUTS.nsave_5d = 1;
OUTPUTS.write_Na00    = '.true.';
OUTPUTS.write_moments = '.true.';
OUTPUTS.write_phi     = '.true.';
OUTPUTS.write_non_lin = '.true.';
OUTPUTS.write_doubleprecision = '.true.';
OUTPUTS.resfile0      = '''results''';
%% Grid parameters
GRID.pmaxe = 1;
GRID.jmaxe = 1;
GRID.pmaxi = 1;
GRID.jmaxi = 1;
GRID.nkr   = 2;
GRID.krmin =-1.0;
GRID.krmax = 1.0;
GRID.nkz   = 2;
GRID.kzmin =-1.0;
GRID.kzmax = 1.0;
%% Model parameters
MODEL.CO      = -2;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
MODEL.NON_LIN = 1;   % Non linear term
MODEL.nu      = 0.001; % collisionality nu*L_perp/Cs0
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
BASIC.nrun                = 1;
BASIC.dt                  = 0.1;
BASIC.tmax                = 0.1;
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
%% Write input file and run HeLaZ
INPUT = write_fort90(OUTPUTS,GRID,MODEL,INITIAL,TIME_INTEGRATION,BASIC);
nproc = 1;
MAKE  = 'cd ..; make; cd wk';
system(MAKE)
EXEC  = ' ../bin/helaz ';
RUN   = ['mpirun -np ' num2str(nproc)];
CMD   = [RUN, EXEC, INPUT];
system(CMD); % RUN HeLaZ

%% Load results
filename = 'results_00.h5';
[Nepj, pi, ji, kr, kz, time] = load_5D_data(filename,'moments_e');
[Nipj, pe, je, kr, kz, time] = load_5D_data(filename,'moments_i');
[Sepj, pe, je, kr, kz, time] = load_5D_data(filename,'Sepj');
[Sipj, pi, ji, kr, kz, time] = load_5D_data(filename,'Sipj');
[phi, kr, kz, time] = load_2D_data(filename,'phi');

%% Compute non linear term with a Matlab method
it = 1; p = 0; j = 0;
for p = 0:1
    for j = 0:1
Sepj_mat = compute_Sapj(p, j, kr, kz, Nepj(:,:,:,:,end), 'e', phi(:,:,end), MODEL, GRID);
    end
end
Sepj_for = squeeze(Sepj(p+1,j+1,:,:,it+1));
% Comparison between Matlab and Fortran method
disp(mean(mean(Sepj_mat - Sepj_for)))
disp(mean(mean(Sepj_for)))
disp(mean(mean(Sepj_mat)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis and basic figures
default_plots_options
SAVEFIG = 1;
filename = 'results_00.h5';
default_plots_options
TITLE  = [TITLE,', $Sapj=',num2str(GRID.kzmin),'$'];
fig = figure;
subplot(121)
plot(kz,Sepj_mat)
hold on
plot(kz,Sepj_for,'--')
subplot(122)
plot(time,squeeze(Sepj(p+1,j+1,1,1,:)))
FIGNAME = 'out';
SIMID = BASIC.SIMID; save_figure