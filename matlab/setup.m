%% _______________________________________________________________________
SIMDIR = ['../results/',SIMID,'/'];
% Grid parameters
GRID.pmaxe = PMAXE;  % Electron Hermite moments
GRID.jmaxe = JMAXE;  % Electron Laguerre moments
GRID.pmaxi = PMAXI;  % Ion Hermite moments
GRID.jmaxi = JMAXI;  % Ion Laguerre moments
GRID.Nx    = NX; % x grid resolution
GRID.Lx    = LX; % x length
GRID.Nexc  = NEXC; % to extend Lx when s>0
GRID.Ny    = NY; % y ''
GRID.Ly    = LY; % y ''
GRID.Nz    = NZ;    % z resolution
GRID.Npol  = NPOL;    % z resolution
if SG; GRID.SG = '.true.'; else; GRID.SG = '.false.';end;
% Geometry
GEOM.geom  = ['''',GEOMETRY,''''];
GEOM.q0    = Q0;    % q factor
GEOM.shear = SHEAR; % shear
GEOM.eps   = EPS;   % inverse aspect ratio
% Model parameters
MODEL.CLOS    = CLOS;
MODEL.NL_CLOS = NL_CLOS;
MODEL.LINEARITY = ['''',LINEARITY,''''];
MODEL.KIN_E   = KIN_E;
if KIN_E; MODEL.KIN_E = '.true.'; else; MODEL.KIN_E = '.false.';end;
MODEL.beta    = BETA;
MODEL.mu_x    = MU_X;
MODEL.mu_y    = MU_Y;
MODEL.N_HD    = N_HD;
MODEL.mu_z    = MU_Z;
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
MODEL.K_E     = 0;        % Electric
MODEL.GradB   = GRADB;      % Magnetic gradient
MODEL.CurvB   = CURVB;      % Magnetic curvature
MODEL.lambdaD = LAMBDAD;
% Collision parameters
COLL.collision_model = ['''',CO,''''];
if (GKCO); COLL.gyrokin_CO = '.true.'; else; COLL.gyrokin_CO = '.false.';end;
if (ABCO); COLL.interspecies = '.true.'; else; COLL.interspecies = '.false.';end;
COLL.mat_file   = '''null''';
switch CO
    case 'SG'
        COLL.mat_file = '''../../../iCa/gk_sugama_P_20_J_10_N_150_kpm_8.0.h5''';
%         COLL.mat_file = '''../../../iCa/gk.hacked_sugama_P_10_J_5_N_150_kpm_8.0.h5''';
%         COLL.mat_file = '''../../../iCa/gk.hacked_sugama_P_4_J_2_N_75_kpm_5.0.h5''';
    case 'LR'
        COLL.mat_file = '''../../../iCa/gk_pitchangle_8_P_20_J_10_N_150_kpm_8.0.h5''';
    case 'LD'
%         COLL.mat_file = '''../../../iCa/gk_coulomb_NFLR_12_P_4_J_2_N_75_kpm_6.0.h5''';
%         COLL.mat_file = '''../../../iCa/gk_coulomb_NFLR_12_P_4_J_2_N_50_kpm_4.0.h5''';
        COLL.mat_file = '''../../../iCa/LDGK_P10_J5_dk_5e-2_km_5_NFLR_12_k2trunc.h5''';
%         COLL.mat_file = '''../../../iCa/LDGK_P10_J5_dk_5e-2_km_5_NFLR_4.h5''';
end
% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = '''RK4''';
INITIAL.INIT_OPT = ['''',INIT_OPT,''''];
INITIAL.ACT_ON_MODES   = ['''',ACT_ON_MODES,''''];
INITIAL.init_background  = BCKGD0;
INITIAL.init_noiselvl    = NOISE0;
INITIAL.iseed            = 42;

% Naming and creating input file
CONAME = CO;
if GKCO
    CONAME = [CONAME,'GK'];
else
    CONAME = [CONAME,'DK'];
end
if ~ABCO
    CONAME = [CONAME,'aa'];
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
coll_ = ['_nu_%1.1e_',CONAME];
coll_   = sprintf(coll_,NU);
% nonlinear
lin_ = [];
if ~LINEARITY; lin_ = '_lin'; end
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
switch GEOMETRY
    case 'circular'
       geo_   = [geo_,'_circ_'];
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
INPUT = write_fort90(OUTPUTS,GRID,GEOM,MODEL,COLL,INITIAL,TIME_INTEGRATION,BASIC);
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
