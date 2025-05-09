%% _______________________________________________________________________
SIMDIR = ['../results/',SIMID,'/'];
% Grid parameters
GRID.pmax = PMAX; % Hermite moments
GRID.jmax = JMAX; % Laguerre moments
GRID.Nx   = NX;   % x grid resolution
GRID.Lx   = LX;   % x length
GRID.Nexc = NEXC; % to extend Lx when s>0
GRID.Ny   = NY;   % y ''
GRID.Ly   = LY;   % y ''
GRID.Nz   = NZ;   % z resolution
if SG; GRID.SG = '.true.'; else; GRID.SG = '.false.';end;
% Geometry
GEOM.geom  = ['''',GEOMETRY,''''];
GEOM.q0    = Q0;    % q factor
GEOM.shear = SHEAR; % shear
GEOM.eps   = EPS;   % inverse aspect ratio
GEOM.kappa   = KAPPA; % elongation
GEOM.s_kappa = S_KAPPA; 
GEOM.delta   = DELTA; % triangularity
GEOM.s_delta = S_DELTA; 
GEOM.zeta    = ZETA;  % squareness
GEOM.s_zeta  = S_ZETA; 
GEOM.parallel_bc  = ['''',PARALLEL_BC,''''];
GEOM.shift_y  = SHIFT_Y;
GEOM.Npol  = NPOL;
if PB_PHASE; GEOM.PB_PHASE = '.true.'; else; GEOM.PB_PHASE = '.false.';end;
% Model parameters
MODEL.LINEARITY = ['''',LINEARITY,''''];
try
    RM_LD_T_EQ;
catch
    RM_LD_T_EQ = 0;
end
if RM_LD_T_EQ; MODEL.RM_LD_T_EQ = '.true.'; else; MODEL.RM_LD_T_EQ = '.false.'; end;
MODEL.Na        = NA;
if ADIAB_E; MODEL.ADIAB_E = '.true.'; else; MODEL.ADIAB_E = '.false.';end;
if ADIAB_I; MODEL.ADIAB_I = '.true.'; else; MODEL.ADIAB_I = '.false.';end;
if MHD_PD;  MODEL.MHD_PD  = '.true.'; else; MODEL.MHD_PD  = '.false.';end;
MODEL.beta    = BETA;
MODEL.ExBrate = EXBRATE;
MODEL.mu_x    = MU_X;
MODEL.mu_y    = MU_Y;
MODEL.N_HD    = N_HD;
try 
    MODEL.N_HDz = N_HDz;
catch
    MODEL.N_HDz = 4;
end
MODEL.mu_z    = MU_Z;
MODEL.HYP_V   = ['''',HYP_V,''''];
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
% gradients L_perp/L_x
MODEL.K_Ni    = K_Ni;       
MODEL.K_Ne    = K_Ne;
MODEL.K_Ti    = K_Ti;    
MODEL.K_Te    = K_Te;    
MODEL.K_gB   = K_gB;      % artificial magnetic gradient  tuner
MODEL.K_cB   = K_cB;      % artificial magnetic curvature tuner
try
    K_mB; K_tB; K_ldB;
catch
    K_mB=1; K_tB=1; K_ldB=1;
end
MODEL.K_mB   = K_mB;      % artificial mirror force   tuner
MODEL.K_tB   = K_tB;      % artificial trapping term  tuner
MODEL.K_ldB  = K_ldB;     % artificial Landau damping tuner
% CLOSURE parameters
CLOSURE.hierarchy_closure = ['''',HRCY_CLOS,''''];
CLOSURE.nonlinear_closure = ['''',NLIN_CLOS,''''];
CLOSURE.dmax              = DMAX;
CLOSURE.nmax              = NMAX;
% Collision parameters
COLL.collision_model = ['''',CO,''''];
if (GKCO); COLL.GK_CO = '.true.'; else; COLL.GK_CO = '.false.';end;
if (ABCO); COLL.INTERSPECIES = '.true.'; else; COLL.INTERSPECIES = '.false.';end;
COLL.mat_file   = 'null';
switch CO
    case 'SG'
        COLL.mat_file = 'gk_sugama_P_20_J_10_N_150_kpm_8.0.h5';
%         COLL.mat_file = 'gk.hacked_sugama_P_10_J_5_N_150_kpm_8.0.h5';
%         COLL.mat_file = 'gk.hacked_sugama_P_4_J_2_N_75_kpm_5.0.h5';
    case 'LR'
        COLL.mat_file = 'gk_pitchangle_8_P_20_J_10_N_150_kpm_8.0.h5';
    case 'LD'
        COLL.mat_file = 'gk_landau_P10_J5_dk_5e-2_km_2.0_NFLR_12.h5';
        % COLL.mat_file = 'gk_landauii_P16_J9_dk_5e-2_km_2.0_NFLR_8.h5';
        % COLL.mat_file = 'gk_landau_P11_J7_dk_5e-2_km_2.0_NFLR_16.h5';
        % COLL.mat_file = 'gk_coulomb_NFLR_12_P_4_J_2_N_50_kpm_4.0.h5';
%         COLL.mat_file = 'LDGKii_P10_J5_dk_5e-2_km_5_NFLR_12_k2trunc.h5';
%         COLL.mat_file = 'LDGKii_P10_J5_dk_5e-2_km_5_NFLR_30.h5';        
%         COLL.mat_file = 'LDGK_P6_J3_dk_5e-2_km_2.5_NFLR_20.h5';        
end
COLL.mat_file = [gyacomodir,'iCa/',COLL.mat_file];
COLL.mat_file = ['''',COLL.mat_file,''''];
COLL.coll_kcut = COLL_KCUT;
% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = ['''',NUMERICAL_SCHEME,''''];
INITIAL.INIT_OPT = ['''',INIT_OPT,''''];
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
CLOSNAME   = [HRCY_CLOS,' dmax=',num2str(DMAX)];
NLCLOSNAME = [NLIN_CLOS,' nmax=',num2str(NMAX)];
% Hermite-Laguerre degrees naming
HLdeg_   = ['_',num2str(PMAX+1),'x',num2str(JMAX+1)];
% temp. dens. drives
drives_ = [];
if abs(K_Ni) > 0; drives_ = [drives_,'_kN_',num2str(K_Ni)]; end;
if abs(K_Ti) > 0; drives_ = [drives_,'_kT_',num2str(K_Ti)]; end;
% collision
coll_ = ['_nu_%1.1e_',CONAME];
coll_   = sprintf(coll_,NU);
% nonlinear
lin_ = [];
if ~LINEARITY; lin_ = '_lin'; end
adiabe_ = [];
if ADIAB_E; adiabe_ = '_adiabe'; end
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
       geo_   = [geo_,'_cir_'];
    case 's-alpha'
       geo_   = [geo_,'_s-a_'];
    case 'miller'
       geo_   = [geo_,'_mil_'];
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
OUTPUTS.dtsave_0d = DTSAVE0D;
OUTPUTS.dtsave_1d = -1;
OUTPUTS.dtsave_2d = DTSAVE2D;
OUTPUTS.dtsave_3d = DTSAVE3D;
OUTPUTS.dtsave_5d = DTSAVE5D;
if W_DOUBLE; OUTPUTS.write_doubleprecision = '.true.'; else; OUTPUTS.write_doubleprecision = '.false.';end;
if W_GAMMA;  OUTPUTS.write_gamma = '.true.'; else; OUTPUTS.write_gamma = '.false.';end;
if W_HF;     OUTPUTS.write_hf    = '.true.'; else; OUTPUTS.write_hf    = '.false.';end;
if W_PHI;    OUTPUTS.write_phi   = '.true.'; else; OUTPUTS.write_phi   = '.false.';end;
if W_NA00;   OUTPUTS.write_Na00  = '.true.'; else; OUTPUTS.write_Na00  = '.false.';end;
if W_NAPJ;   OUTPUTS.write_Napj  = '.true.'; else; OUTPUTS.write_Napj  = '.false.';end;
if W_SAPJ;   OUTPUTS.write_Sapj  = '.true.'; else; OUTPUTS.write_Sapj  = '.false.';end;
if W_DENS;   OUTPUTS.write_dens  = '.true.'; else; OUTPUTS.write_dens  = '.false.';end;
if W_FVEL;   OUTPUTS.write_fvel  = '.true.'; else; OUTPUTS.write_fvel  = '.false.';end;
if W_TEMP;   OUTPUTS.write_temp  = '.true.'; else; OUTPUTS.write_temp  = '.false.';end;
OUTPUTS.job2load    = JOB2LOAD;
%% Create directories
if ~exist(SIMDIR, 'dir')
   mkdir(SIMDIR)
end
if ~exist(BASIC.RESDIR, 'dir')
mkdir(BASIC.RESDIR)
end
% if ~exist(BASIC.MISCDIR, 'dir')
% mkdir(BASIC.MISCDIR)
% end
%% Compile and WRITE input file
INPUT = write_fort90(OUTPUTS,GRID,GEOM,MODEL,CLOSURE,COLL,INITIAL,TIME_INTEGRATION,BASIC);
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
