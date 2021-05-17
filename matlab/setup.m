%% ________________________________________________________________________
SIMDIR = ['../results/',SIMID,'/'];
% Grid parameters
GRID.pmaxe = PMAXE;  % Electron Hermite moments
GRID.jmaxe = JMAXE;  % Electron Laguerre moments
GRID.pmaxi = PMAXI;  % Ion Hermite moments
GRID.jmaxi = JMAXI;  % Ion Laguerre moments
GRID.Nr    = N; % r grid resolution
GRID.Lr    = L; % r length
GRID.Nz    = N * (1-KREQ0) + KREQ0; % z ''
GRID.Lz    = L * (1-KREQ0); % z ''
GRID.kpar  = KPAR;

% Model parameters
MODEL.CO      = CO;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
MODEL.CLOS    = CLOS;
MODEL.NL_CLOS = NL_CLOS;
if NON_LIN; MODEL.NON_LIN = '.true.'; else; MODEL.NON_LIN = '.false.';end;
MODEL.mu      = MU;
MODEL.mu_p    = MU_P;
MODEL.mu_j    = MU_J;
MODEL.nu      = NU; % hyper diffusive coefficient nu for HW
% temperature ratio T_a/T_e
MODEL.tau_e   = TAU;
MODEL.tau_i   = TAU;
% mass ratio sqrt(m_a/m_i)
MODEL.sigma_e = 0.0233380;
MODEL.sigma_i = 1.0;
% charge q_a/e
MODEL.q_e     =-1.0;
MODEL.q_i     = 1.0;
if MODEL.q_e == 0; SIMID = [SIMID,'_i']; end;
% gradients L_perp/L_x
MODEL.eta_n   = ETAN;        % source term kappa for HW
MODEL.eta_T   = ETAT;        % Temperature
MODEL.eta_B   = ETAB;        % Magnetic
MODEL.lambdaD = LAMBDAD;
% if A0KH ~= 0; SIMID = [SIMID,'_Nz_',num2str(L/2/pi*KR0KH),'_A_',num2str(A0KH)]; end;
% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = '''RK4''';
if INIT_PHI; INITIAL.init_noisy_phi = '.true.'; else; INITIAL.init_noisy_phi = '.false.';end;
INITIAL.init_background  = 0.0e-5;
INITIAL.init_noiselvl = NOISE0;
INITIAL.iseed             = 42;
INITIAL.selfmat_file = '''null''';
INITIAL.eimat_file = '''null''';
INITIAL.iemat_file = '''null''';
if (CO == -3) % Write matrice filename for Full Coulomb DK
    cmat_pmaxe = 25;
    cmat_jmaxe = 18;
    cmat_pmaxi = 25;
    cmat_jmaxi = 18;
    INITIAL.selfmat_file = ...
        ['''../../../iCa/self_Coll_GKE_0_GKI_0_ESELF_1_ISELF_1_Pmaxe_',num2str(cmat_pmaxe),...
        '_Jmaxe_',num2str(cmat_jmaxe),'_Pmaxi_',num2str(cmat_pmaxi),'_Jmaxi_',...
        num2str(cmat_jmaxi),'_pamaxx_10.h5'''];
    INITIAL.eimat_file = ...
        ['''../../../iCa/ei_Coll_GKE_0_GKI_0_ETEST_1_EBACK_1_Pmaxe_',num2str(cmat_pmaxe),...
        '_Jmaxe_',num2str(cmat_jmaxe),'_Pmaxi_',num2str(cmat_pmaxi),'_Jmaxi_',...
        num2str(cmat_jmaxi),'_pamaxx_10_tau_1.0000_mu_0.0233.h5'''];
    INITIAL.iemat_file = ...
        ['''../../../iCa/ie_Coll_GKE_0_GKI_0_ITEST_1_IBACK_1_Pmaxe_',num2str(cmat_pmaxe),...
        '_Jmaxe_',num2str(cmat_jmaxe),'_Pmaxi_',num2str(cmat_pmaxi),'_Jmaxi_',...
        num2str(cmat_jmaxi),'_pamaxx_10_tau_1.0000_mu_0.0233.h5'''];
elseif (CO == -2) % Write matrice filename for DK Sugama
    cmat_pmaxe = 10;
    cmat_jmaxe = 5;
    cmat_pmaxi = 10;
    cmat_jmaxi = 5;
    INITIAL.selfmat_file = ...
        ['''../../../iCa/self_Coll_GKE_0_GKI_0_ESELF_3_ISELF_3_Pmaxe_',num2str(cmat_pmaxe),...
        '_Jmaxe_',num2str(cmat_jmaxe),'_Pmaxi_',num2str(cmat_pmaxi),'_Jmaxi_',...
        num2str(cmat_jmaxi),'_JE_12.h5'''];
    INITIAL.eimat_file = ...
        ['''../../../iCa/ei_Coll_GKE_0_GKI_0_ETEST_3_EBACK_3_Pmaxe_',num2str(cmat_pmaxe),...
        '_Jmaxe_',num2str(cmat_jmaxe),'_Pmaxi_',num2str(cmat_pmaxi),'_Jmaxi_',...
        num2str(cmat_jmaxi),'_JE_12_tau_1.0000_mu_0.0233.h5'''];
    INITIAL.iemat_file = ...
        ['''../../../iCa/ie_Coll_GKE_0_GKI_0_ITEST_3_IBACK_3_Pmaxe_',num2str(cmat_pmaxe),...
        '_Jmaxe_',num2str(cmat_jmaxe),'_Pmaxi_',num2str(cmat_pmaxi),'_Jmaxi_',...
        num2str(cmat_jmaxi),'_JE_12_tau_1.0000_mu_0.0233.h5'''];
elseif (CO == 2) % Write matrice filename for Sugama GK
    cmat_pmaxe = 10;
    cmat_jmaxe = 5;
    cmat_pmaxi = 10;
    cmat_jmaxi = 5;
    INITIAL.selfmat_file = ...
        ['''../../../iCa/self_Coll_GKE_1_GKI_1_ESELF_3_ISELF_3_Pmaxe_',num2str(cmat_pmaxe),...
        '_Jmaxe_',num2str(cmat_jmaxe),'_Pmaxi_',num2str(cmat_pmaxi),'_Jmaxi_',...
        num2str(cmat_jmaxi),'_JE_12_'''];
    INITIAL.eimat_file = ...
        ['''../../../iCa/ei_Coll_GKE_1_GKI_1_ETEST_3_EBACK_3_Pmaxe_',num2str(cmat_pmaxe),...
        '_Jmaxe_',num2str(cmat_jmaxe),'_Pmaxi_',num2str(cmat_pmaxi),'_Jmaxi_',...
        num2str(cmat_jmaxi),'_JE_12_tau_1.0000_mu_0.0233_'''];
    INITIAL.iemat_file = ...
        ['''../../../iCa/ie_Coll_GKE_1_GKI_1_ITEST_3_IBACK_3_Pmaxe_',num2str(cmat_pmaxe),...
        '_Jmaxe_',num2str(cmat_jmaxe),'_Pmaxi_',num2str(cmat_pmaxi),'_Jmaxi_',...
        num2str(cmat_jmaxi),'_JE_12_tau_1.0000_mu_0.0233_'''];
elseif (CO == 3) % Full Coulomb GK
    disp('Warning, FCGK not implemented yet')
elseif (CO == -1) % DGDK
    disp('Warning, DGDK not debugged')
end

% Naming and creating input file
if    (CO == -3); CONAME = 'FCDK';
elseif(CO == -2); CONAME = 'SGDK';
elseif(CO == -1); CONAME = 'DGDK';
elseif(CO ==  0); CONAME = 'LB';
elseif(CO ==  1); CONAME = 'DGGK';
elseif(CO ==  2); CONAME = 'SGGK';
elseif(CO ==  3); CONAME = 'FCGK';
end
if    (CLOS == 0); CLOSNAME = 'Trunc.';
elseif(CLOS == 1); CLOSNAME = 'Clos. 1';
elseif(CLOS == 2); CLOSNAME = 'Clos. 2';
end
if (PMAXE == PMAXI) && (JMAXE == JMAXI)
    degngrad   = ['P_',num2str(PMAXE),'_J_',num2str(JMAXE)];
else
    degngrad   = ['Pe_',num2str(PMAXE),'_Je_',num2str(JMAXE),...
        '_Pi_',num2str(PMAXI),'_Ji_',num2str(JMAXI)];
end
degngrad = [degngrad,'_eta_',num2str(ETAB/ETAN),'_nu_%0.0e_',...
        CONAME,'_CLOS_',num2str(CLOS),'_mu_%0.0e'];
    
degngrad   = sprintf(degngrad,[NU,MU]);
if ~NON_LIN; degngrad = ['lin_',degngrad]; end
resolution = [num2str(GRID.Nr),'x',num2str(GRID.Nz/2),'_'];
gridname   = ['L_',num2str(L),'_'];
if (exist('PREFIX','var') == 0); PREFIX = []; end;
if (exist('SUFFIX','var') == 0); SUFFIX = []; end;
PARAMS = [PREFIX,resolution,gridname,degngrad,SUFFIX];
BASIC.RESDIR = [SIMDIR,PARAMS,'/'];
BASIC.PARAMS = PARAMS;
BASIC.SIMID  = SIMID;
BASIC.nrun       = 1e8;
BASIC.dt         = DT;
BASIC.tmax       = TMAX;    %time normalized to 1/omega_pe
BASIC.maxruntime = str2num(CLUSTER.TIME(1:2))*3600 ...
                   + str2num(CLUSTER.TIME(4:5))*60 ...
                   + str2num(CLUSTER.TIME(7:8));
% Outputs parameters
if RESTART; BASIC.RESTART = '.true.'; else; BASIC.RESTART = '.false.';end;
OUTPUTS.nsave_0d = floor(1.0/SPS0D/DT);
OUTPUTS.nsave_1d = -1;
OUTPUTS.nsave_2d = floor(1.0/SPS2D/DT);
OUTPUTS.nsave_5d = floor(1.0/SPS5D/DT);
OUTPUTS.nsave_cp = floor(1.0/SPSCP/DT);
if W_DOUBLE; OUTPUTS.write_doubleprecision = '.true.'; else; OUTPUTS.write_doubleprecision = '.false.';end;
if W_GAMMA;  OUTPUTS.write_gamma = '.true.'; else; OUTPUTS.write_gamma = '.false.';end;
if W_PHI;    OUTPUTS.write_phi   = '.true.'; else; OUTPUTS.write_phi   = '.false.';end;
if W_NA00;   OUTPUTS.write_Na00  = '.true.'; else; OUTPUTS.write_Na00  = '.false.';end;
if W_NAPJ;   OUTPUTS.write_Napj  = '.true.'; else; OUTPUTS.write_Napj  = '.false.';end;
if W_SAPJ;   OUTPUTS.write_Sapj  = '.true.'; else; OUTPUTS.write_Sapj  = '.false.';end;
OUTPUTS.resfile0    = '''outputs''';
OUTPUTS.rstfile0    = '''checkpoint''';
OUTPUTS.job2load    = JOB2LOAD;
%% Create directories
if ~exist(SIMDIR, 'dir')
   mkdir(SIMDIR)
end
if ~exist(BASIC.RESDIR, 'dir')
mkdir(BASIC.RESDIR)
end
%% Compile and WRITE input file
INPUT = write_fort90(OUTPUTS,GRID,MODEL,INITIAL,TIME_INTEGRATION,BASIC);
nproc = 1;
MAKE  = 'cd ..; make; cd wk';
system(MAKE);
%%
disp(['Set up ',SIMID]);
disp([resolution,gridname,degngrad]);
if RESTART
	disp(['- restarting from JOBNUM = ',num2str(JOB2LOAD)]); else
	disp(['- starting from T = 0']);
end
