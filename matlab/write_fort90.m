function [INPUT] = write_fort90(OUTPUTS,GRID,MODEL,INITIAL,TIME_INTEGRATION,BASIC)
% Write the input script "fort.90" with desired parameters
INPUT = 'fort.90';
fid = fopen(INPUT,'wt');

fprintf(fid,'&BASIC\n');
fprintf(fid,['  nrun   = ', num2str(BASIC.nrun),'\n']);
fprintf(fid,['  dt     = ', num2str(BASIC.dt),'\n']);
fprintf(fid,['  tmax   = ', num2str(BASIC.tmax),'\n']);
fprintf(fid,['  RESTART    = ', BASIC.RESTART,'\n']);
fprintf(fid,['  maxruntime = ', num2str(BASIC.maxruntime),'\n']);
fprintf(fid,'/\n');

fprintf(fid,'&GRID\n');
fprintf(fid,['  pmaxe  =', num2str(GRID.pmaxe),'\n']);
fprintf(fid,['  jmaxe  = ', num2str(GRID.jmaxe),'\n']);
fprintf(fid,['  pmaxi  = ', num2str(GRID.pmaxi),'\n']);
fprintf(fid,['  jmaxi  = ', num2str(GRID.jmaxi),'\n']);
fprintf(fid,['  Nx     = ', num2str(GRID.Nx),'\n']);
fprintf(fid,['  Lx     = ', num2str(GRID.Lx),'\n']);
fprintf(fid,['  Ny     = ', num2str(GRID.Ny),'\n']);
fprintf(fid,['  Ly     = ', num2str(GRID.Ly),'\n']);
fprintf(fid,['  Nz     = ', num2str(GRID.Nz),'\n']);
fprintf(fid,['  q0     = ', num2str(GRID.q0),'\n']);
fprintf(fid,['  shear  = ', num2str(GRID.shear),'\n']);
fprintf(fid,['  eps    = ', num2str(GRID.eps),'\n']);
fprintf(fid,'/\n');

fprintf(fid,'&OUTPUT_PAR\n');
fprintf(fid,['  nsave_0d = ', num2str(OUTPUTS.nsave_0d),'\n']);
fprintf(fid,['  nsave_1d = ', num2str(OUTPUTS.nsave_1d),'\n']);
fprintf(fid,['  nsave_2d = ', num2str(OUTPUTS.nsave_2d),'\n']);
fprintf(fid,['  nsave_3d = ', num2str(OUTPUTS.nsave_3d),'\n']);
fprintf(fid,['  nsave_5d = ', num2str(OUTPUTS.nsave_5d),'\n']);
fprintf(fid,['  nsave_cp = ', num2str(OUTPUTS.nsave_cp),'\n']);
fprintf(fid,['  write_doubleprecision = ', OUTPUTS.write_doubleprecision,'\n']);
fprintf(fid,['  write_gamma = ', OUTPUTS.write_gamma,'\n']);
fprintf(fid,['  write_phi   = ', OUTPUTS.write_phi,'\n']);
fprintf(fid,['  write_Na00  = ', OUTPUTS.write_Na00,'\n']);
fprintf(fid,['  write_Napj  = ', OUTPUTS.write_Napj,'\n']);
fprintf(fid,['  write_Sapj  = ', OUTPUTS.write_Sapj,'\n']);
fprintf(fid,['  write_dens  = ', OUTPUTS.write_dens,'\n']);
fprintf(fid,['  write_temp  = ', OUTPUTS.write_temp,'\n']);
fprintf(fid,['  resfile0      = ', OUTPUTS.resfile0,'\n']);
fprintf(fid,['  rstfile0      = ', OUTPUTS.rstfile0,'\n']);
fprintf(fid,['  job2load      = ', num2str(OUTPUTS.job2load),'\n']);
fprintf(fid,'/\n');

fprintf(fid,'&MODEL_PAR\n');
fprintf(fid,'  ! Collisionality\n');
fprintf(fid,['  CO      = ', num2str(MODEL.CO),'\n']);
fprintf(fid,['  CLOS    = ', num2str(MODEL.CLOS),'\n']);
fprintf(fid,['  NL_CLOS = ', num2str(MODEL.NL_CLOS),'\n']);
fprintf(fid,['  NON_LIN = ', MODEL.NON_LIN,'\n']);
fprintf(fid,['  mu      = ', num2str(MODEL.mu),'\n']);
fprintf(fid,['  mu_p    = ', num2str(MODEL.mu_p),'\n']);
fprintf(fid,['  mu_j    = ', num2str(MODEL.mu_j),'\n']);
fprintf(fid,['  nu      = ', num2str(MODEL.nu),'\n']);
fprintf(fid,['  tau_e   = ', num2str(MODEL.tau_e),'\n']);
fprintf(fid,['  tau_i   = ', num2str(MODEL.tau_i),'\n']);
fprintf(fid,['  sigma_e = ', num2str(MODEL.sigma_e),'\n']);
fprintf(fid,['  sigma_i = ', num2str(MODEL.sigma_i),'\n']);
fprintf(fid,['  q_e     = ', num2str(MODEL.q_e),'\n']);
fprintf(fid,['  q_i     = ', num2str(MODEL.q_i),'\n']);
fprintf(fid,['  eta_n   = ', num2str(MODEL.eta_n),'\n']);
fprintf(fid,['  eta_T   = ', num2str(MODEL.eta_T),'\n']);
fprintf(fid,['  eta_B   = ', num2str(MODEL.eta_B),'\n']);
fprintf(fid,['  lambdaD = ', num2str(MODEL.lambdaD),'\n']);
fprintf(fid,'/\n');

fprintf(fid,'&INITIAL_CON\n');
fprintf(fid,['  INIT_NOISY_PHI    = ', INITIAL.init_noisy_phi,'\n']);
fprintf(fid,['  INIT_ZF       = ', num2str(INITIAL.INIT_ZF),'\n']);
fprintf(fid,['  WIPE_TURB     = ', INITIAL.wipe_turb,'\n']);
fprintf(fid,['  INIT_BLOB     = ', INITIAL.init_blob,'\n']);
fprintf(fid,['  init_background  = ', num2str(INITIAL.init_background),'\n']);
fprintf(fid,['  init_noiselvl = ', num2str(INITIAL.init_noiselvl),'\n']);
fprintf(fid,['  iseed         = ', num2str(INITIAL.iseed),'\n']);
fprintf(fid,['  mat_file      = ', INITIAL.mat_file,'\n']);
fprintf(fid,'/\n');

fprintf(fid,'&TIME_INTEGRATION_PAR\n');
fprintf(fid,['  numerical_scheme = ', TIME_INTEGRATION.numerical_scheme,'\n']);
fprintf(fid,'/');

fclose(fid);
system(['cp fort.90 ',BASIC.RESDIR,'/.']);
end
