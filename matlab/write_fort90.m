function [INPUT] = write_fort90(OUTPUTS,GRID,GEOM,MODEL,COLL,INITIAL,TIME_INTEGRATION,BASIC)
% Write the input script "fort.90" with desired parameters
INPUT = ['fort_',sprintf('%2.2d',OUTPUTS.job2load+1),'.90'];
fid = fopen(INPUT,'wt');

fprintf(fid,'&BASIC\n');
fprintf(fid,['  nrun   = ', num2str(BASIC.nrun),'\n']);
fprintf(fid,['  dt     = ', num2str(BASIC.dt),'\n']);
fprintf(fid,['  tmax   = ', num2str(BASIC.tmax),'\n']);
fprintf(fid,['  maxruntime = ', num2str(BASIC.maxruntime),'\n']);
fprintf(fid,'/\n');

fprintf(fid,'&GRID\n');
fprintf(fid,['  pmaxe  = ', num2str(GRID.pmaxe),'\n']);
fprintf(fid,['  jmaxe  = ', num2str(GRID.jmaxe),'\n']);
fprintf(fid,['  pmaxi  = ', num2str(GRID.pmaxi),'\n']);
fprintf(fid,['  jmaxi  = ', num2str(GRID.jmaxi),'\n']);
fprintf(fid,['  Nx     = ', num2str(GRID.Nx),'\n']);
fprintf(fid,['  Lx     = ', num2str(GRID.Lx),'\n']);
fprintf(fid,['  Ny     = ', num2str(GRID.Ny),'\n']);
fprintf(fid,['  Ly     = ', num2str(GRID.Ly),'\n']);
fprintf(fid,['  Nz     = ', num2str(GRID.Nz),'\n']);
fprintf(fid,['  Npol   = ', num2str(GRID.Npol),'\n']);
fprintf(fid,['  SG     = ',           GRID.SG,'\n']);
fprintf(fid,['  Nexc   = ', num2str(GRID.Nexc),'\n']);
fprintf(fid,'/\n');

fprintf(fid,'&GEOMETRY\n');
fprintf(fid,['  geom   = ', GEOM.geom,'\n']);
fprintf(fid,['  q0     = ', num2str(GEOM.q0),'\n']);
fprintf(fid,['  shear  = ', num2str(GEOM.shear),'\n']);
fprintf(fid,['  eps    = ', num2str(GEOM.eps),'\n']);
fprintf(fid,['  kappa  = ', num2str(GEOM.kappa),'\n']);
fprintf(fid,['  delta  = ', num2str(GEOM.delta),'\n']);
fprintf(fid,['  zeta   = ', num2str(GEOM.zeta),'\n']);
fprintf(fid,['  parallel_bc = ', GEOM.parallel_bc,'\n']);
fprintf(fid,'/\n');

fprintf(fid,'&OUTPUT_PAR\n');
fprintf(fid,['  nsave_0d = ', num2str(OUTPUTS.nsave_0d),'\n']);
fprintf(fid,['  nsave_1d = ', num2str(OUTPUTS.nsave_1d),'\n']);
fprintf(fid,['  nsave_2d = ', num2str(OUTPUTS.nsave_2d),'\n']);
fprintf(fid,['  nsave_3d = ', num2str(OUTPUTS.nsave_3d),'\n']);
fprintf(fid,['  nsave_5d = ', num2str(OUTPUTS.nsave_5d),'\n']);
fprintf(fid,['  write_doubleprecision = ', OUTPUTS.write_doubleprecision,'\n']);
fprintf(fid,['  write_gamma = ', OUTPUTS.write_gamma,'\n']);
fprintf(fid,['  write_hf    = ', OUTPUTS.write_hf,'\n']);
fprintf(fid,['  write_phi   = ', OUTPUTS.write_phi,'\n']);
fprintf(fid,['  write_Na00  = ', OUTPUTS.write_Na00,'\n']);
fprintf(fid,['  write_Napj  = ', OUTPUTS.write_Napj,'\n']);
fprintf(fid,['  write_Sapj  = ', OUTPUTS.write_Sapj,'\n']);
fprintf(fid,['  write_dens  = ', OUTPUTS.write_dens,'\n']);
fprintf(fid,['  write_temp  = ', OUTPUTS.write_temp,'\n']);
fprintf(fid,['  job2load    = ', num2str(OUTPUTS.job2load),'\n']);
fprintf(fid,'/\n');

fprintf(fid,'&MODEL_PAR\n');
fprintf(fid,'  ! Collisionality\n');
fprintf(fid,['  CLOS    = ', num2str(MODEL.CLOS),'\n']);
fprintf(fid,['  NL_CLOS = ', num2str(MODEL.NL_CLOS),'\n']);
fprintf(fid,['  LINEARITY = ', MODEL.LINEARITY,'\n']);
fprintf(fid,['  KIN_E   = ', MODEL.KIN_E,'\n']);
fprintf(fid,['  mu_x    = ', num2str(MODEL.mu_x),'\n']);
fprintf(fid,['  mu_y    = ', num2str(MODEL.mu_y),'\n']);
fprintf(fid,['  N_HD    = ', num2str(MODEL.N_HD),'\n']);
fprintf(fid,['  mu_z    = ', num2str(MODEL.mu_z),'\n']);
fprintf(fid,['  HYP_V   = ', MODEL.HYP_V,'\n']);
fprintf(fid,['  mu_p    = ', num2str(MODEL.mu_p),'\n']);
fprintf(fid,['  mu_j    = ', num2str(MODEL.mu_j),'\n']);
fprintf(fid,['  nu      = ', num2str(MODEL.nu),'\n']);
fprintf(fid,['  tau_e   = ', num2str(MODEL.tau_e),'\n']);
fprintf(fid,['  tau_i   = ', num2str(MODEL.tau_i),'\n']);
fprintf(fid,['  sigma_e = ', num2str(MODEL.sigma_e),'\n']);
fprintf(fid,['  sigma_i = ', num2str(MODEL.sigma_i),'\n']);
fprintf(fid,['  q_e     = ', num2str(MODEL.q_e),'\n']);
fprintf(fid,['  q_i     = ', num2str(MODEL.q_i),'\n']);
fprintf(fid,['  K_Ne    = ', num2str(MODEL.K_Ne),'\n']);
fprintf(fid,['  K_Te    = ', num2str(MODEL.K_Te),'\n']);
fprintf(fid,['  K_Ni    = ', num2str(MODEL.K_Ni),'\n']);
fprintf(fid,['  K_Ti    = ', num2str(MODEL.K_Ti),'\n']);
fprintf(fid,['  k_gB   = ', num2str(MODEL.k_gB),'\n']);
fprintf(fid,['  k_cB   = ', num2str(MODEL.k_cB),'\n']);
fprintf(fid,['  lambdaD = ', num2str(MODEL.lambdaD),'\n']);
fprintf(fid,['  beta    = ', num2str(MODEL.beta),'\n']);
fprintf(fid,'/\n');

fprintf(fid,'&COLLISION_PAR\n');
fprintf(fid,['  collision_model = ', COLL.collision_model,'\n']);
fprintf(fid,['  gyrokin_CO      = ', COLL.gyrokin_CO,'\n']);
fprintf(fid,['  interspecies    = ', COLL.interspecies,'\n']);
fprintf(fid,['  mat_file        = ', COLL.mat_file,'\n']);
fprintf(fid,['  collision_kcut  = ', num2str(COLL.coll_kcut),'\n']);
fprintf(fid,'/\n');


fprintf(fid,'&INITIAL_CON\n');
fprintf(fid,['  INIT_OPT      = ', INITIAL.INIT_OPT,'\n']);
fprintf(fid,['  ACT_ON_MODES  = ', INITIAL.ACT_ON_MODES,'\n']);
fprintf(fid,['  init_background  = ', num2str(INITIAL.init_background),'\n']);
fprintf(fid,['  init_noiselvl = ', num2str(INITIAL.init_noiselvl),'\n']);
fprintf(fid,['  iseed         = ', num2str(INITIAL.iseed),'\n']);
fprintf(fid,'/\n');

fprintf(fid,'&TIME_INTEGRATION_PAR\n');
fprintf(fid,['  numerical_scheme = ', TIME_INTEGRATION.numerical_scheme,'\n']);
fprintf(fid,'/');

fclose(fid);
system(['cp fort*.90 ',BASIC.RESDIR,'/.']);
end
