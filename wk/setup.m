%% ________________________________________________________________________
% Grid parameters
GRID.pmaxe = PMAXE;  % Electron Hermite moments
GRID.jmaxe = JMAXE;  % Electron Laguerre moments 
GRID.pmaxi = PMAXI;  % Ion Hermite moments
GRID.jmaxi = JMAXI;  % Ion Laguerre moments
GRID.Nr    = N * (1-KREQ0) + KREQ0; % r grid resolution
GRID.Lr    = L * (1-KREQ0); % r length
GRID.Nz    = N; % z ''
GRID.Lz    = L; % z ''
GRID.kpar  = KPAR;
% Model parameters
MODEL.CO      = CO;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
if 0;      MODEL.DK      = '.true.'; else; MODEL.DK      = '.false.';end;
if NON_LIN; MODEL.NON_LIN = '.true.'; else; MODEL.NON_LIN = '.false.';end;
MODEL.mu      = MU;
MODEL.nu      = NU; % hyper diffusive coefficient nu for HW
% temperature ratio T_a/T_e
MODEL.tau_e   = TAU;
MODEL.tau_i   = TAU;
% mass ratio sqrt(m_a/m_i)
MODEL.sigma_e = 0.0233380;
MODEL.sigma_i = 1.0;
% charge q_a/e
MODEL.q_e     =-1.0 * (1-NO_E);
MODEL.q_i     = 1.0;
if MODEL.q_e == 0; SIMID = [SIMID,'_i']; end;
% gradients L_perp/L_x
MODEL.eta_n   = ETAN;        % source term kappa for HW
MODEL.eta_T   = ETAT;        % Temperature
MODEL.eta_B   = ETAB;        % Magnetic
MODEL.lambdaD = LAMBDAD;
% background phi drive for Kelvin-Helmholtz instability
MODEL.kr0KH   = KR0KH;
MODEL.A0KH    = A0KH;
if A0KH ~= 0; SIMID = [SIMID,'_Nz_',num2str(L/2/pi*KR0KH),'_A_',num2str(A0KH)]; end;
% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = '''RK4''';
INITIAL.only_Na00         = '.false.';
INITIAL.initback_moments  = 0.0e-5;
INITIAL.initnoise_moments = NOISE0;
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
% Naming and creating input file
if    (CO == 0); CONAME = '';
elseif(CO == -1); CONAME = 'FC';
elseif(CO == -2); CONAME = 'DG';
end
params   = ['Pe_',num2str(PMAXE),'_Je_',num2str(JMAXE),...
    '_Pi_',num2str(PMAXI),'_Ji_',num2str(JMAXI),...
    '_nB_',num2str(ETAB),'_nN_',num2str(ETAN),'_nu_%0.0e_',...
    CONAME,'_mu_%0.0e_'];
params   = sprintf(params,[NU,MU]);
if ~NON_LIN; params = ['lin_',params]; end
resolution = ['_',num2str(GRID.Nr),'x',num2str(GRID.Nz/2),'_'];
gridname   = ['L_',num2str(L),'_'];
BASIC.SIMID = [SIMID,resolution,gridname,params];
BASIC.nrun       = 1e8;
BASIC.dt         = DT;   
BASIC.tmax       = TMAX;    %time normalized to 1/omega_pe
% Outputs parameters
if RESTART; BASIC.RESTART = '.true.'; else; BASIC.RESTART = '.false.';end;
OUTPUTS.nsave_0d = 0;
OUTPUTS.nsave_1d = 0;
OUTPUTS.nsave_2d = floor(1.0/SPS2D/DT);
OUTPUTS.nsave_5d = floor(1.0/SPS5D/DT);
OUTPUTS.write_Na00    = '.true.';
OUTPUTS.write_moments = '.true.';
OUTPUTS.write_phi     = '.true.';
OUTPUTS.write_non_lin = OUTPUTS.write_moments;
if NON_LIN == 0; OUTPUTS.write_non_lin = '.false.'; end;
OUTPUTS.write_doubleprecision = '.true.';
OUTPUTS.resfile0      = ['''',BASIC.SIMID,''''];
OUTPUTS.rstfile0      = ['''','../checkpoint/cp_',BASIC.SIMID,''''];
OUTPUTS.job2load      = JOB2LOAD;
%% Compile and write input file
INPUT = write_fort90(OUTPUTS,GRID,MODEL,INITIAL,TIME_INTEGRATION,BASIC);
nproc = 1;
MAKE  = 'cd ..; make; cd wk';
system(MAKE);
%%
disp(['Set up ', BASIC.SIMID]);
if RESTART
	disp(['- restarting from JOBNUM = ',num2str(JOB2LOAD)]); else
	disp(['- starting from T = 0']);
end
