%% ________________________________________________________________________
SIMDIR = ['../results/',SIMID,'/'];
% Grid parameters
GRID.pmaxe = PMAXE;  % Electron Hermite moments
GRID.jmaxe = JMAXE;  % Electron Laguerre moments
GRID.pmaxi = PMAXI;  % Ion Hermite moments
GRID.jmaxi = JMAXI;  % Ion Laguerre moments
GRID.Nx    = NX; % x grid resolution
GRID.Lx    = LX; % x length
GRID.Ny    = NY; % y ''
GRID.Ly    = LY; % y ''
GRID.Nz    = NZ;    % z resolution
GRID.q0    = Q0;    % q factor
GRID.shear = SHEAR; % shear
GRID.eps   = EPS;   % inverse aspect ratio

% Model parameters
MODEL.CO      = CO;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
MODEL.CLOS    = CLOS;
MODEL.NL_CLOS = NL_CLOS;
if NON_LIN; MODEL.NON_LIN = '.true.'; else; MODEL.NON_LIN = '.false.';end;
MODEL.KIN_E = KIN_E;
if KIN_E; MODEL.KIN_E = '.true.'; else; MODEL.KIN_E = '.false.';end;
MODEL.mu      = MU;
MODEL.mu_p    = MU_P;
MODEL.mu_j    = MU_J;
MODEL.nu      = NU; % hyper diffusive coefficient nu for HW
% temperature ratio T_a/T_e
MODEL.tau_e   = TAU;
MODEL.tau_i   = TAU;
% mass ratio sqrt(m_a/m_i)
MODEL.sigma_e = SIGMA_E;
MODEL.sigma_i = 1.0;
% charge q_a/e
MODEL.q_e     =-1.0;
MODEL.q_i     = 1.0;
if MODEL.q_e == 0; SIMID = [SIMID,'_i']; end;
% gradients L_perp/L_x
MODEL.K_n     = K_N;        % source term kappa for HW
MODEL.K_T     = K_T;        % Temperature
MODEL.K_E     = K_E;        % Electric
MODEL.GradB   = GRADB;      % Magnetic gradient
MODEL.CurvB   = CURVB;      % Magnetic curvature
MODEL.lambdaD = LAMBDAD;
% if A0KH ~= 0; SIMID = [SIMID,'_Nz_',num2str(L/2/pi*KX0KH),'_A_',num2str(A0KH)]; end;
% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = '''RK4''';
if (INIT_PHI); INITIAL.init_noisy_phi = '.true.'; else; INITIAL.init_noisy_phi = '.false.';end;
INITIAL.INIT_ZF = INIT_ZF;
INITIAL.wipe_turb = WIPE_TURB;
INITIAL.wipe_zf   = WIPE_ZF;
if (INIT_BLOB); INITIAL.init_blob = '.true.'; else; INITIAL.init_blob = '.false.';end;
INITIAL.init_background  = (INIT_ZF>0)*ZF_AMP + BCKGD0;
INITIAL.init_noiselvl = NOISE0;
INITIAL.iseed             = 42;
INITIAL.mat_file   = '''null''';
if (abs(CO) == 2) %Sugama operator
    INITIAL.mat_file = ['''../../../iCa/gk_sugama_P_20_J_10_N_150_kpm_8.0.h5'''];
elseif (abs(CO) == 3) %pitch angle operator
    INITIAL.mat_file = ['''../../../iCa/gk_pitchangle_8_P_20_J_10_N_150_kpm_8.0.h5'''];
elseif (CO == 4) % Full Coulomb GK
%     INITIAL.mat_file = ['''../../../iCa/gk_coulomb_NFLR_12_P_4_J_2_N_50_kpm_4.0.h5'''];
    INITIAL.mat_file = ['''../../../iCa/gk_coulomb_NFLR_12_P_4_J_2_N_75_kpm_6.0.h5'''];
%     INITIAL.mat_file = ['''../../../iCa/gk_coulomb_NFLR_6_P_4_J_2_N_50_kpm_4.0.h5'''];
%     INITIAL.mat_file = ['''../../../iCa/gk_coulomb_NFLR_6_P_4_J_2_N_75_kpm_6.0.h5'''];
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
% Hermite-Laguerre degrees naming
if (PMAXE == PMAXI) && (JMAXE == JMAXI)
    HLdeg_   = ['_',num2str(PMAXE+1),'x',num2str(JMAXE+1)];
else
    HLdeg_   = ['_Pe_',num2str(PMAXE+1),'_Je_',num2str(JMAXE+1),...
        '_Pi_',num2str(PMAXI+1),'_Ji_',num2str(JMAXI+1)];
end
% temp. dens. drives
drives_ = [];
if abs(K_N) > 0; drives_ = [drives_,'_kN_',num2str(K_N)]; end;
if abs(K_T) > 0; drives_ = [drives_,'_kT_',num2str(K_T)]; end;
% collision
coll_ = ['_nu_%0.0e_',CONAME];
coll_   = sprintf(coll_,NU);
% nonlinear
lin_ = [];
if ~NON_LIN; lin_ = '_lin'; end
adiabe_ = [];
if ~KIN_E; adiabe_ = '_adiabe'; end
% resolution and boxsize
res_ = [num2str(GRID.Nx),'x',num2str(GRID.Ny)];
if  (LX ~= LY)
    geo_   = ['_Lx_',num2str(LX),'_Ly_',num2str(LY)];
else
    geo_   = ['_L_',num2str(LX)];
end
if (NZ > 1)  %3D case
    res_ = [res_,'x',num2str(NZ)];
    if abs(Q0) > 0
        geo_   = [geo_,'_q0_',num2str(Q0)];
    end
    if abs(EPS) > 0
       geo_   = [geo_,'_e_',num2str(EPS)];
    end
    if abs(SHEAR) > 0
       geo_   = [geo_,'_s_',num2str(SHEAR)];
    end
end
% put everything together in the param character chain
u_ = '_'; % underscore variable
PARAMS = [res_,HLdeg_,geo_,drives_,coll_,lin_,adiabe_];
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
OUTPUTS.nsave_0d = floor(1.0/SPS0D/DT);
OUTPUTS.nsave_1d = -1;
OUTPUTS.nsave_2d = floor(1.0/SPS2D/DT);
OUTPUTS.nsave_3d = floor(1.0/SPS3D/DT);
OUTPUTS.nsave_5d = floor(1.0/SPS5D/DT);
if W_DOUBLE; OUTPUTS.write_doubleprecision = '.true.'; else; OUTPUTS.write_doubleprecision = '.false.';end;
if W_GAMMA;  OUTPUTS.write_gamma = '.true.'; else; OUTPUTS.write_gamma = '.false.';end;
if W_HF;     OUTPUTS.write_hf    = '.true.'; else; OUTPUTS.write_hf    = '.false.';end;
if W_PHI;    OUTPUTS.write_phi   = '.true.'; else; OUTPUTS.write_phi   = '.false.';end;
if W_NA00;   OUTPUTS.write_Na00  = '.true.'; else; OUTPUTS.write_Na00  = '.false.';end;
if W_NAPJ;   OUTPUTS.write_Napj  = '.true.'; else; OUTPUTS.write_Napj  = '.false.';end;
if W_SAPJ;   OUTPUTS.write_Sapj  = '.true.'; else; OUTPUTS.write_Sapj  = '.false.';end;
if W_DENS;   OUTPUTS.write_dens  = '.true.'; else; OUTPUTS.write_dens  = '.false.';end;
if W_TEMP;   OUTPUTS.write_temp  = '.true.'; else; OUTPUTS.write_temp  = '.false.';end;
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
% system(MAKE);
%%
disp(['Set up ',SIMID]);
disp([res_,geo_,HLdeg_]);
if JOB2LOAD>=0
	disp(['- restarting from JOBNUM = ',num2str(JOB2LOAD)]); else
	disp(['- starting from T = 0']);
end
