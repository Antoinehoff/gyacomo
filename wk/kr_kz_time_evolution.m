%% Run linear simulation on a kr,kz grid and compute the non linear term afterwards
clear all; close all;
BASIC.SIMID = 'kr_kz_time_evolution'; % Name of the simulations
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
GRID.pmaxe = 8;
GRID.jmaxe = 4;
GRID.pmaxi = 8;
GRID.jmaxi = 4;
GRID.nkr   = 16;
GRID.krmin = 0.1;
GRID.krmax = 1.;
GRID.nkz   = 16;
GRID.kzmin = 0.1;
GRID.kzmax = 1.0;
%% Model parameters
MODEL.CO      = -2;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
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
MODEL.eta_n   = 0.0;        % Density
MODEL.eta_T   = 0.0;        % Temperature
MODEL.eta_B   = 1.0;        % Magnetic
% Debye length
MODEL.lambdaD = 0.0;
%% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = '''RK4''';
BASIC.nrun                = 100000;
BASIC.dt                  = 0.01;
BASIC.tmax                = 100.0;
INITIAL.initback_moments  = 0.01;
INITIAL.initnoise_moments = 0.0;
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
EXEC  = ' ../bin/helaz ';
RUN   = ['mpirun -np ' num2str(nproc)];
CMD   = [RUN, EXEC, INPUT];
system(CMD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis and basic figures
default_plots_options
SAVEFIG = 1;
filename = 'results_00.h5';
default_plots_options
TITLE  = [TITLE,', $k_z=',num2str(GRID.kzmin),'$'];

%% Load HeLaZ moments

moment = 'Ni00';

kr     = h5read(filename,['/data/var2d/' moment '/coordkr']);
kz     = h5read(filename,['/data/var2d/' moment '/coordkz']);
time   = h5read(filename,'/data/var2d/time');
Ni00   = zeros(numel(kr),numel(kz),numel(time));
for it = 1:numel(time)
    tmp          = h5read(filename,['/data/var2d/', moment,'/', num2str(it,'%06d')]);
    Ni00(:,:,it) = tmp.real + 1i * tmp.imaginary;
end

%% Gif creation
disp('Creating gif..');
skip    = 50;    % To skip some frames
delay   = 0.05; % Speed of the movie (smaller for faster)
t       = time(1:skip:end);
N       = squeeze(Ni00(:,:,1:skip:end));
H       = zeros(GRID.nkr, GRID.nkz, numel(t));
h       = zeros(GRID.nkr, GRID.nkz, numel(t));
Pad = 1.0; % Padding
for it = 1:numel(t)
    H(:,:,it) = abs(N(:,:,it));
    F         = abs(ifft2(H,floor(Pad*GRID.nkr),floor(Pad*GRID.nkz)));
    h(:,:,it) = F(1:GRID.nkr, 1:GRID.nkr);
end
GIFNAME = 'Ni00';
create_gif(kr, kz, t, H, BASIC, GRID, MODEL, delay, GIFNAME)
% GIFNAME = 'Mi00';
% create_gif(kr, kz, t, h, BASIC, GRID, MODEL, delay, GIFNAME)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
