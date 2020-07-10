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
OUTPUTS.write_Ni00    = '.true.';
OUTPUTS.write_moments = '.true.';
OUTPUTS.write_phi     = '.true.';
OUTPUTS.write_doubleprecision = '.true.';
OUTPUTS.resfile0      = '''results''';
%% Grid parameters
GRID.pmaxe = 2;
GRID.jmaxe = 2;
GRID.pmaxi = 2;
GRID.jmaxi = 2;
GRID.nkr   = 10;
GRID.krmin = 0.0;
GRID.krmax = 0.;
GRID.nkz   = 10;
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
MODEL.eta_n   = 1.0;        % Density
MODEL.eta_T   = 0.0;        % Temperature
MODEL.eta_B   = 0.5;        % Magnetic
% Debye length
MODEL.lambdaD = 0.0;
%% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = '''RK4''';
BASIC.nrun                = 100000;
BASIC.dt                  = 0.05;
BASIC.tmax                = 1.0;
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
EXEC  = ' ../bin/helaz ';
RUN   = ['mpirun -np ' num2str(nproc)];
CMD   = [RUN, EXEC, INPUT];
system(CMD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run MOLI for each grid point and save the moments evolution
kr_scan = linspace(GRID.krmin,GRID.krmax,GRID.nkr);
kz_scan = linspace(GRID.kzmin,GRID.kzmax,GRID.nkz);
Nmtot   = (GRID.pmaxe+1) * (GRID.jmaxe+1) + (GRID.pmaxi+1) * (GRID.jmaxi+1);
Nt      = floor(BASIC.tmax/BASIC.dt)+1;
Napj_2D = zeros(GRID.nkr,GRID.nkz,Nt,Nmtot);
params.RK4 = 1;
for ikr = 1:GRID.nkr
    for ikz = 1:GRID.nkz
        kr = kr_scan(ikr); % used by MOLI_time_solver_2D directly
        kz = kz_scan(ikz); 
        run ../matlab/MOLI_time_solver_2D
        Napj_2D(ikr,ikz,:,:) = results.Napj;
    end
end
%% Compute non linear term

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis and basic figures
default_plots_options
SAVEFIG = 1;
filename = 'results_00.h5';
default_plots_options
TITLE  = [TITLE,', $k_z=',num2str(GRID.kzmin),'$'];
if params.RK4; MOLIsolvername = 'RK4';
else;          MOLIsolvername = 'ode15i';
end
bare = @(p_,j_) (GRID.jmaxe+1)*p_ + j_ + 1;
bari = @(p_,j_) bare(GRID.pmaxe, GRID.jmaxe) + (GRID.jmaxi+1)*p_ + j_ + 1;

%% Convert MOLI moments into HeLaZ format
Nepj_MOLI = zeros(GRID.pmaxe+1, GRID.jmaxe+1,GRID.nkr,GRID.nkz,Nt);
for ip = 1:GRID.pmaxe+1
    for ij = 1:GRID.jmaxe+1
        Nepj_MOLI(ip,ij,:,:,:) = Napj_2D(:,:,:,bare(ip-1,ij-1));
    end
end
Nipj_MOLI = zeros(GRID.pmaxi+1, GRID.jmaxi+1,GRID.nkr,GRID.nkz,Nt);
for ip = 1:GRID.pmaxi+1
    for ij = 1:GRID.jmaxi+1
        Nipj_MOLI(ip,ij,:,:,:) = Napj_2D(:,:,:,bari(ip-1,ij-1));
    end
end
%% Load HeLaZ moments

moment = 'moments_i';

kr     = h5read(filename,['/data/var5d/' moment '/coordkr']);
kz     = h5read(filename,['/data/var5d/' moment '/coordkz']);
time   = h5read(filename,'/data/var5d/time');
Nipj_HeLaZ   = zeros(GRID.pmaxi+1, GRID.jmaxi+1,numel(kr),numel(kz),numel(time));
for it = 1:numel(time)
    tmp          = h5read(filename,['/data/var5d/', moment,'/', num2str(it,'%06d')]);
    Nipj_HeLaZ(:,:,:,:,it) = tmp.real + 1i * tmp.imaginary;
end

moment = 'moments_e';

kr     = h5read(filename,['/data/var5d/' moment '/coordkr']);
kz     = h5read(filename,['/data/var5d/' moment '/coordkz']);
time   = h5read(filename,'/data/var5d/time');
Nepj_HeLaZ   = zeros(GRID.pmaxe+1, GRID.jmaxe+1,numel(kr),numel(kz),numel(time));
for it = 1:numel(time)
    tmp          = h5read(filename,['/data/var5d/', moment,'/', num2str(it,'%06d')]);
    Nepj_HeLaZ(:,:,:,:,it) = tmp.real + 1i * tmp.imaginary;
end

%% Check coherence of linear results
disp('MOLI - HeLaZ differences :');
disp(sum(sum(sum(sum(sum(Nipj_MOLI-Nipj_HeLaZ)))))...
     / sum(sum(sum(sum(sum(Nipj_MOLI))))));
 
%% Non linear term computation
%get phi
phiHeLaZ  = zeros(numel(timephi),numel(kr),numel(kz));
for it = 1:numel(time)
    tmp  = h5read(filename,['/data/var2d/phi/' num2str(it,'%06d')]);
    phiHeLaZ(it,:,:) = tmp.real + 1i * tmp.imaginary;
end

compute_Sapj

