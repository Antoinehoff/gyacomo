%% ________________________________________________________________________
SIMDIR = ['../results/',SIMID,'/'];
% Grid parameters
GRID.pmaxe = PMAXE;  % Electron Hermite moments
GRID.jmaxe = JMAXE;  % Electron Laguerre moments
GRID.pmaxi = PMAXI;  % Ion Hermite moments
GRID.jmaxi = JMAXI;  % Ion Laguerre moments
GRID.Nx    = N; % x grid resolution
GRID.Lx    = L; % x length
GRID.Ny    = N * (1-KXEQ0) + KXEQ0; % y ''
GRID.Ly    = L * (1-KXEQ0); % y ''
GRID.Nz    = Nz;    % z resolution
GRID.q0    = q0;    % q factor
GRID.shear = shear; % shear
GRID.eps   = eps;   % inverse aspect ratio

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
% if A0KH ~= 0; SIMID = [SIMID,'_Nz_',num2str(L/2/pi*KX0KH),'_A_',num2str(A0KH)]; end;
% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = '''RK4''';
if (INIT_PHI); INITIAL.init_noisy_phi = '.true.'; else; INITIAL.init_noisy_phi = '.false.';end;
INITIAL.INIT_ZF = INIT_ZF;
if (WIPE_TURB); INITIAL.wipe_turb = '.true.'; else; INITIAL.wipe_turb = '.false.';end;
if (INIT_BLOB); INITIAL.init_blob = '.true.'; else; INITIAL.init_blob = '.false.';end;
INITIAL.init_background  = (INIT_ZF>0)*ZF_AMP;
INITIAL.init_noiselvl = NOISE0;
INITIAL.iseed             = 42;
INITIAL.mat_file   = '''null''';
if (abs(CO) == 2) %Sugama operator
    INITIAL.mat_file = ['''../../../iCa/gk_sugama_P_20_J_10_N_150_kpm_8.0.h5'''];
elseif (abs(CO) == 3) %pitch angle operator
    INITIAL.mat_file = ['''../../../iCa/gk_pitchangle_8_P_20_J_10_N_150_kpm_8.0.h5'''];
elseif (CO == 4) % Full Coulomb GK
    INITIAL.mat_file = ['''../../../iCa/gk_coulomb_P_6_J_3_N_150_kpm_8.0_NFLR_4.h5'''];
%     INITIAL.mat_file = ['''../../../iCa/gk_coulomb_P_6_J_3_N_150_kpm_8.0_NFLR_k2.h5'''];
elseif (CO == -1) % DGDK
    disp('Warning, DGDK not debugged')
end

% Naming and creating input file
switch abs(CO)
    case 0; CONAME = 'LB';
    case 1; CONAME = 'DG';
    case 2; CONAME = 'SG';
    case 3; CONAME = 'PA';
    case 4; CONAME = 'FC';
    otherwise; CONAME ='UK';
end
if (CO <= 0); CONAME = [CONAME,'DK'];
else;         CONAME = [CONAME,'GK'];
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
        CONAME,'_mu_%0.0e'];

degngrad   = sprintf(degngrad,[NU,MU]);
if ~NON_LIN; degngrad = ['lin_',degngrad]; end
if (Nz == 1)
    resolution = [num2str(GRID.Nx),'x',num2str(GRID.Ny/2),'_'];
    gridname   = ['L_',num2str(L),'_'];
else
    resolution = [num2str(GRID.Nx),'x',num2str(GRID.Ny/2),'x',num2str(GRID.Nz),'_'];
    gridname   = ['L_',num2str(L),'_q0_',num2str(q0),'_'];
end
if (exist('PREFIX','var') == 0); PREFIX = []; end;
if (exist('SUFFIX','var') == 0); SUFFIX = []; end;
PARAMS = [PREFIX,resolution,gridname,degngrad,SUFFIX];
BASIC.RESDIR  = [SIMDIR,PARAMS,'/'];
BASIC.MISCDIR = ['/misc/HeLaZ_outputs/',SIMDIR(4:end),PARAMS,'/'];
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
OUTPUTS.nsave_3d = floor(1.0/SPS3D/DT);
OUTPUTS.nsave_5d = floor(1.0/SPS5D/DT);
OUTPUTS.nsave_cp = floor(1.0/SPSCP/DT);
if W_DOUBLE; OUTPUTS.write_doubleprecision = '.true.'; else; OUTPUTS.write_doubleprecision = '.false.';end;
if W_GAMMA;  OUTPUTS.write_gamma = '.true.'; else; OUTPUTS.write_gamma = '.false.';end;
if W_PHI;    OUTPUTS.write_phi   = '.true.'; else; OUTPUTS.write_phi   = '.false.';end;
if W_NA00;   OUTPUTS.write_Na00  = '.true.'; else; OUTPUTS.write_Na00  = '.false.';end;
if W_NAPJ;   OUTPUTS.write_Napj  = '.true.'; else; OUTPUTS.write_Napj  = '.false.';end;
if W_SAPJ;   OUTPUTS.write_Sapj  = '.true.'; else; OUTPUTS.write_Sapj  = '.false.';end;
if W_DENS;   OUTPUTS.write_dens  = '.true.'; else; OUTPUTS.write_dens  = '.false.';end;
if W_TEMP;   OUTPUTS.write_temp  = '.true.'; else; OUTPUTS.write_temp  = '.false.';end;
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
if ~exist(BASIC.MISCDIR, 'dir')
mkdir(BASIC.MISCDIR)
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
