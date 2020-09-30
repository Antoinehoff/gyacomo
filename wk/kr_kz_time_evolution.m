% Run linear/nonlin simulation on a kr,kz grid and create gifs on the
% fourier and real space
clear all; close all;
BASIC.SIMID = 'kr_kz_time_evolution'; % Name of the simulations
addpath(genpath('../matlab')) % ... add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% outputs options
OUTPUTS.nsave_0d = 0;
OUTPUTS.nsave_1d = 0;
OUTPUTS.nsave_2d = 10;
OUTPUTS.nsave_5d = 0;
OUTPUTS.write_Ni00    = '.true.';
OUTPUTS.write_moments = '.true.';
OUTPUTS.write_phi     = '.true.';
OUTPUTS.write_non_lin = '.false.';
OUTPUTS.write_doubleprecision = '.true.';
OUTPUTS.resfile0      = '''results''';
OUTPUTS.rstfile0      = '''restart''';
%% Grid parameters
GRID.pmaxe = 08;  % Electron Hermite moments
GRID.jmaxe = 04;  % Electron Laguerre moments
GRID.pmaxi = 08;  % Ion Hermite moments
GRID.jmaxi = 04;  % Ion Laguerre moments
GRID.nkr   = 128; % kr grid resolution
GRID.krmin =-2.0; % kr minimum value
GRID.krmax = 2.0; % kr maximal value
GRID.nkz   = 128; % kz ''
GRID.kzmin =-2.0; %    ''
GRID.kzmax = 2.0; %    ''
GRID.Pad   = 2.0; % Zero padding for dealiasing (Mx = Pad * nkx)
%% Model parameters
MODEL.CO      = -2;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
MODEL.NON_LIN = '.false.';   % Non linear term
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
BASIC.dt                  = 0.01;
BASIC.tmax                = 200.0;    %time normalized to 1/omega_pe
INITIAL.RESTART           = '.true.';
INITIAL.backup_file      = '''restart''';
INITIAL.only_Na00         = '.false.';
INITIAL.initback_moments  = 1.0e-4;
INITIAL.initnoise_moments = 5.0e-5;
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
system(MAKE);
EXEC  = ' ../bin/helaz ';
RUN   = ['mpirun -np ' num2str(nproc)];
CMD   = [RUN, EXEC, INPUT];
%system(CMD);
%% load results
filename = [OUTPUTS.resfile0(2:end-1),'_00.h5'];
[Ni00, kr, kz, time] = load_2D_data(filename, 'Ni00');
[ phi, kr, kz, time] = load_2D_data(filename, 'phi');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis and basic figures
default_plots_options
SAVEFIG = 1;
default_plots_options
%% Gif creation
if strcmp(MODEL.NON_LIN,'.true.'); LINEARITY = 'nl';
else; LINEARITY = 'lin'; end;
disp('Creating gif..');
disp('');
skip    = 5;    % To skip some frames
delay   = 0.03; % Speed of the movie (smaller for faster)
t       = time(1:skip:end);
N       = squeeze(Ni00(:,:,1:skip:end));
H       = zeros(GRID.nkr, GRID.nkz, numel(t));
h       = zeros(GRID.nkr, GRID.nkz, numel(t));
Pad = 2.0; % Padding
for it = 1:numel(t)
    H(:,:,it) = abs(N(:,:,it));
    F         = real(ifft2(N(:,:,it),ceil(Pad*GRID.nkr)+1,ceil(Pad*GRID.nkz)+1));
    h(:,:,it) = F(1:GRID.nkr, 1:GRID.nkr);
end
%%
GIFNAME = [num2str(GRID.nkr),'x',num2str(GRID.nkz),'-',LINEARITY,'-Ni00'];
create_gif(kr, kz, t, log(H), BASIC, GRID, MODEL, delay, GIFNAME, false)

%%
Nr = ceil(GRID.nkr/2); Nz = ceil(GRID.nkz/2);
h_ = h;
h_(h_==0) = nan;
GIFNAME = [num2str(GRID.nkr),'x',num2str(GRID.nkz),'-',LINEARITY,'-Mi00'];
create_gif(kr(Nr:end), kz(Nr:end), t, h_(Nr:end,Nz:end,:), BASIC, GRID, MODEL, delay, GIFNAME, false)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
