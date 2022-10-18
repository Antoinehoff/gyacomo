%% Directory of the simulation
gyacomodir = pwd; gyacomodir = gyacomodir(1:end-2); % get code directory
% if 1% Local results
resdir ='';
resdir ='';
resdir ='';
resdir ='';
% resdir ='debug/ppj_init';
%% nu = 5e-1
% Sugama
% resdir ='Hallenbert_nu_5e-01/200x32_5x3_L_120_kN_1.5_kT_0.375_nu_5e-01_SGGK';% also in 7x4
% resdir ='Hallenbert_nu_5e-01/200x32_5x3_L_120_kN_1.6_kT_0.4_nu_5e-01_SGGK';
% resdir ='Hallenbert_nu_5e-01/200x32_7x4_L_120_kN_1.7_kT_0.425_nu_5e-01_SGGK';% also in 7x4
% resdir ='Hallenbert_nu_5e-01/200x32_5x3_L_120_kN_1.8_kT_0.45_nu_5e-01_SGGK';
% resdir ='Hallenbert_nu_5e-01/200x32_5x3_L_120_kN_1.9_kT_0.475_nu_5e-01_SGGK';%also in 7x4

%% nu = 1e-1
% Landau
% resdir ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.5_kT_0.375_nu_1e-01_LDGK';
% resdir ='Hallenbert_nu_1e-01/150x50_5x3_L_120_kN_1.6_kT_0.4_nu_1e-01_LDGK';
% resdir ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.6_kT_0.4_nu_1e-01_LDGK';
% resdir ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.7_kT_0.425_nu_1e-01_LDGK';
% resdir ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.8_kT_0.45_nu_1e-01_LDGK';
% resdir ='Hallenbert_nu_1e-01/150x50_5x3_L_120_kN_1.8_kT_0.45_nu_1e-01_LDGK';
% resdir ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.9_kT_0.475_nu_1e-01_LDGK';

% Sugama
% resdir ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.5_kT_0.375_nu_1e-01_SGGK';
% resdir ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.7_kT_0.425_nu_1e-01_SGGK';
% resdir ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.8_kT_0.45_nu_1e-01_SGGK';
% resdir ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.9_kT_0.475_nu_1e-01_SGGK';

% Dougherty
% resdir ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.5_kT_0.375_nu_1e-01_DGGK';
% resdir ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.6_kT_0.4_nu_1e-01_DGGK';
% resdir ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.7_kT_0.425_nu_1e-01_DGGK';
% resdir ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.8_kT_0.45_nu_1e-01_DGGK';

%% nu = 5e-2
% resdir ='Hallenbert_nu_5e-02/200x32_11x6_L_120_kN_1.8_kT_0.45_nu_5e-02_SGGK';%For GENE benchmark % to analyse (added HD)
% resdir ='Hallenbert_nu_5e-02/200x32_11x6_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_SGGK';

% testing various NL closures
% resdir ='Hallenbert_nu_5e-02/200x32_7x4_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_SGGK';

% resdir ='Hallenbert_nu_5e-02/200x32_5x3_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_SGDK';
% resdir ='Hallenbert_nu_5e-02/200x32_11x6_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_SGDK';

% resdir ='Hallenbert_nu_5e-02/256x64_5x3_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_SGDK';
% resdir ='Hallenbert_nu_5e-02/256x64_11x6_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_SGDK';
% resdir ='Hallenbert_nu_5e-02/200x32_21x3_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_SGDK';
% resdir ='Hallenbert_nu_5e-02/200x32_17x9_L_120_kN_1.8_kT_0.45_nu_5e-02_SGDK';

% resdir ='Hallenbert_nu_5e-02/200x32_5x3_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_DGGK';
% resdir ='Hallenbert_nu_5e-02/200x32_11x6_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_DGGK';
% resdir ='Hallenbert_nu_5e-02/128x32_5x3_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_FCGK';

%% nu = 1e-2
% Landau
% resdir ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.5_kT_0.375_nu_1e-02_LDGK';
% resdir ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.6_kT_0.4_nu_1e-02_LDGK';
% resdir ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.7_kT_0.425_nu_1e-02_LDGK';
% resdir ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.8_kT_0.45_nu_1e-02_LDGK';
% resdir ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.9_kT_0.475_nu_1e-02_LDGK';

% Sugama
% resdir ='kobayashi_2015_fig1/150x150_5x3_L_100_kN_1.4_nu_5e-03_SGGK';
% resdir ='Hallenbert_nu_1e-02/200x32_7x4_L_120_kN_1.5_kT_0.375_nu_1e-02_SGGK';
% resdir ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.5_kT_0.375_nu_1e-02_SGGK';
% resdir ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.6_kT_0.4_nu_1e-02_SGGK';
% resdir ='Hallenbert_nu_1e-02/200x32_7x4_L_120_kN_1.6_kT_0.4_nu_1e-02_SGGK';
% resdir ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.7_kT_0.425_nu_1e-02_SGGK';
% resdir ='Hallenbert_nu_1e-02/300x64_5x3_L_120_kN_1.7_kT_0.425_nu_1e-02_SGGK';
% resdir ='Hallenbert_nu_1e-02/200x32_11x6_L_120_kN_1.7_kT_0.425_nu_1e-02_SGGK';
% resdir ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.8_kT_0.45_nu_1e-02_SGGK'; % To analyse (added HD)
% resdir ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.9_kT_0.475_nu_1e-02_SGGK';
% resdir ='Hallenbert_nu_1e-02/200x32_11x6_L_120_kN_1.9_kT_0.475_nu_1e-02_SGGK';
% Dougherty
% resdir ='Hallenbert_nu_1e-02/200x32_7x4_L_120_kN_1.5_kT_0.375_nu_1e-02_DGGK';
% resdir ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.6_kT_0.4_nu_1e-02_DGGK';
% resdir ='Hallenbert_nu_1e-02/200x32_8x5_L_120_kN_1.6_kT_0.4_nu_0e+00_DGGK';
% resdir ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.8_kT_0.45_nu_1e-02_DGGK';

%% nu = 5e-3
% resdir ='Hallenbert_nu_5e-03/200x32_5x3_Lx_120_Ly_60_kN_1.8_eta_0.25_nuSG_5e-03_muxy_5e-2';
% resdir ='Hallenbert_nu_5e-03/200x32_5x3_Lx_120_Ly_60_kN_1.8_eta_0.25_nuSG_5e-03_mux_5e-2_muy_6e-1';
% resdir ='Hallenbert_nu_5e-03/200x32_11x6_Lx_120_Ly_60_kN_1.8_eta_0.25_nuSG_5e-03_mux_5e-2_muy_6e-1';
% resdir ='Hallenbert_nu_5e-03/200x32_11x6_Lx_120_Ly_60_kN_1.8_eta_0.25_nuSG_5e-03_muxy_5e-2';

%% nu = 0

% resdir ='Hallenbert_fig2a/200x32_21x11_Lx_120_Ly_60_kN_1.6_eta_0.4_nu_0_muxy_1e-2';
% resdir ='Hallenbert_fig2a/200x32_11x6_Lx_120_Ly_60_kN_1.6_eta_0.4_nu_0_muxy_1e-2';
% resdir ='Hallenbert_fig2a/200x32_5x3_Lx_120_Ly_60_kN_1.6_eta_0.4_nu_0_muxy_1e-2';
% resdir ='Hallenbert_fig2a/200x32_5x3_Lx_120_Ly_60_kN_1.6_eta_0.4_nuDGGK_0.1_muxy_1e-2';
% resdir ='Hallenbert_fig2a/200x32_11x6_Lx_120_Ly_60_kN_1.6_eta_0.4_nuDGGK_0.1_muxy_1e-2';
% resdir ='Hallenbert_fig2a/200x32_5x3_Lx_120_Ly_60_kN_1.6_eta_0.4_nuSGGK_0.1_muxy_1e-2';
% resdir ='Hallenbert_fig2a/200x32_11x6_Lx_120_Ly_60_kN_1.6_eta_0.4_nuSGGK_0.1_muxy_1e-2';
% resdir ='Hallenbert_fig2a/200x32_5x3_Lx_120_Ly_60_kN_1.6_eta_0.4_nuLDGK_0.1_muxy_1e-2';
% resdir ='Hallenbert_fig2a/200x32_5x3_Lx_120_Ly_60_kN_1.6_eta_0.4_nuLRGK_0.1_muxy_1e-2';


% resdir ='Hallenbert_fig2b/200x32_11x6_Lx_240_Ly_120_kN_2.5_eta_0.25_nu_0_muxy_1e-1';
% resdir ='Hallenbert_fig2b/200x32_5x3_Lx_240_Ly_120_kN_2.5_eta_0.25_nu_0_muxy_1e-1';
% resdir ='Hallenbert_fig2b/200x32_5x3_Lx_240_Ly_120_kN_2.5_eta_0.25_nuSGGK_0.1_muxy_1e-1';
% resdir ='Hallenbert_fig2b/200x32_5x3_Lx_240_Ly_120_kN_2.5_eta_0.25_nuDGGK_0.1_muxy_1e-1';
% resdir ='Hallenbert_fig2b/200x32_5x3_Lx_240_Ly_120_kN_2.5_eta_0.25_nuLDGK_0.1_muxy_1e-1';
% resdir ='Hallenbert_fig2b/200x32_5x3_Lx_240_Ly_120_kN_2.5_eta_0.25_nuLRGK_0.1_muxy_1e-1';

% resdir ='Hallenbert_fig2c/200x32_11x6_Lx_120_Ly_60_kN_2.0_eta_0.25_nu_0_muxy_5e-2';
% resdir ='Hallenbert_fig2c/200x32_5x3_Lx_120_Ly_60_kN_2.0_eta_0.25_nu_0_muxy_5e-2';
% resdir ='Hallenbert_fig2c/200x32_5x3_Lx_120_Ly_60_kN_2.0_eta_0.25_nuSGGK_0.1_muxy_5e-2';
% resdir ='Hallenbert_fig2c/200x32_11x6_Lx_120_Ly_60_kN_2.0_eta_0.25_nuSGGK_0.1_muxy_5e-2';
% resdir ='Hallenbert_fig2c/200x32_5x3_Lx_120_Ly_60_kN_2.0_eta_0.25_nuDGGK_0.1_muxy_5e-2';
% resdir ='Hallenbert_fig2c/200x32_11x6_Lx_120_Ly_60_kN_2.0_eta_0.25_nuDGGK_0.1_muxy_5e-2';
% resdir ='Hallenbert_fig2c/200x32_5x3_Lx_120_Ly_60_kN_2.0_eta_0.25_nuLDGK_0.1_muxy_5e-2';
% resdir ='Hallenbert_fig2c/200x32_5x3_Lx_120_Ly_60_kN_2.0_eta_0.25_nuLRGK_0.1_muxy_5e-2';
% resdir ='Hallenbert_fig2c/200x32_9x5_Lx_120_Ly_60_kN_2.0_eta_0.25_nuLRGK_0.1_muxy_5e-2';

%% Transport scan
% resdir = 'nu_0.1_transport_scan/colless_kn_1.7_to_2.0';
% resdir = 'nu_0.1_transport_scan/colless_kn_2.1_to_2.5';

% resdir = 'nu_0.1_transport_scan/LB_kn_2.0';

% resdir = 'nu_0.1_transport_scan/DG_kn_1.8_to_2.1';
% resdir = 'nu_0.1_transport_scan/DG_kn_2.2_to_2.5';
% resdir = 'nu_0.1_transport_scan/DG_conv_kN_1.9';

% resdir = 'nu_0.1_transport_scan/SG_kn_1.7_to_2.0';
% resdir = 'nu_0.1_transport_scan/SG_10x5_conv_test';
% resdir = 'nu_0.1_transport_scan/SG_kn_2.2_to_2.5';

% resdir = 'nu_0.1_transport_scan/LD_kn_2.0_to_2.5';
% resdir = 'nu_0.1_transport_scan/LD_kn_1.7_to_2.5';

% resdir = 'nu_0.1_transport_scan/LR_kn_1.7_to_2.0';
% resdir = 'nu_0.1_transport_scan/LR_kn_2.1_to_2.5';

% resdir = 'nu_0.1_transport_scan/colless_kn_2.2_Lx1.5';
% resdir = 'nu_0.1_transport_scan/colless_kn_2.2_HD';

% resdir = 'nu_0.1_transport_scan/colless_kn_1.6_HD';

% resdir = 'nu_0.1_transport_scan/large_box_kN_2.1_nu_0.1';
% resdir = 'nu_0.1_transport_scan/large_box_kN_2.0_nu_0.1';

% resdir = 'predator_prey_nu_scan/DG_Kn_1.7_nu_0.01';

% resdir = 'ZF_damping_linear_nu_0_20x10_kn_1.6_GK/LR_4x2_nu_0.1';
% resdir = 'ZF_damping_nu_0_20x10_kn_1.6_GK/HSG_4x2_nu_0.1';
% resdir = 'ZF_damping_nu_0_5x3_kn_2.5_GK/LR_4x2_nu_0.1';
% resdir = 'hacked_sugama/hacked_B_kn_1.6_200x32_L_120x60_nu_0.1';

% resdir = 'shearless_cyclone/200x32x24_5x4_Lx_120_Ly_60_q0_1.4_e_0.18_kN_2.22_kT_6.9_nuLR_0.01_adiab_e';
% resdir = 'shearless_cyclone/no_sg_128x32x36_6x3_Lx_120_Ly_60_q0_1.4_e_0.18_kN_2.22_kT_6.9_adiab_e';
% resdir = 'shearless_cyclone/sgrid_128x64x32_4x2_Lx_100_Ly_120_q0_1.4_e_0.18_kN_2.22_kT_6.96_adiab_e';
% resdir = 'shearless_cyclone/sgrid_128x64x32_4x2_Lx_100_Ly_120_q0_1.4_e_0.18_kN_1.78_kT_5.52_adiab_e';
% resdir = 'linear_shearless_cyclone/4_2_cyclone_1.0';
% resdir = 'linear_shearless_cyclone/test_fmom';
% else% Marconi results
% resdir ='';
% resdir ='';
% resdir ='';fd
% resdir ='';
% resdir ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/simulation_A_new/300x300_5x3_L_120_kN_1.6667_nu_1e-01_SGGK/out.txt';
% % resdir ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/simulation_A/300x150_L_120_P_8_J_4_eta_0.6_nu_1e-01_SGGK_mu_0e+00/out.txt';
% % BASIC.RESDIR      = ['../',resdir(46:end-8),'/'];
% MISCDIR = ['/misc/HeLaZ_outputs/',resdir(46:end-8),'/'];
% end

%% ZPINCH rerun
% resdir ='Zpinch_rerun/Kn_2.5_200x48x5x3';
% resdir ='Zpinch_rerun/Kn_2.5_256x128x5x3';
% resdir ='Zpinch_rerun/Kn_2.5_312x196x5x3_Lx_400_Ly_200';
% resdir ='Zpinch_rerun/Kn_2.5_256x64x5x3';
% resdir ='Zpinch_rerun/Kn_2.0_200x48x9x5_large_box';
% resdir ='Zpinch_rerun/Kn_2.0_256x64x9x5_Lx_240_Ly_120';
% resdir ='Zpinch_rerun/Kn_1.6_256x128x7x4';
% resdir ='Zpinch_rerun/Kn_1.6_200x64x11x6';
% resdir ='Zpinch_rerun/Kn_1.6_200x64x11x6_conv';
% resdir ='Zpinch_rerun/Kn_1.6_200x64x11x6_mu_0.5';
% resdir ='Zpinch_rerun/Kn_1.6_256x128x21x11';

%% nu scan
% resdir = 'Zpinch_rerun/kN_2.2_coll_scan_128x48x5x3';
% resdir = 'Zpinch_rerun/Ultra_HD_312x196x5x3';
% resdir = 'Zpinch_rerun/UHD_nu_001_LDGK';
% resdir = 'Zpinch_rerun/UHD_nu_01_LDGK';
% resdir = 'Zpinch_rerun/UHD_nu_1_LDGK';
% resdir ='Zpinch_rerun/kN_1.7_SGGK_conv_200x32x7x3_nu_0.01';
% resdir ='Zpinch_rerun/kN_1.7_LDGKii_200x32x7x3_nu_scan';
% resdir ='Zpinch_rerun/nu_0.1_LDGKii_200x48x7x4_kN_scan';
% resdir ='Zpinch_rerun/nu_0.1_FCGK_200x48x5x3_kN_scan';
% resdir ='Zpinch_rerun/nu_0.1_SGGK_200x48x5x3_kN_scan';
% resdir = 'Zpinch_rerun/kN_1.7_FCGK_200x32x5x3_nu_scan';
% % resdir = 'Zpinch_rerun/kN_1.7_SGGK_200x32x7x4_nu_scan';
resdir = 'Zpinch_rerun/kN_2.2_SGGK_200x32x5x3_nu_scan';
% resdir = 'Zpinch_rerun/kN_1.7_SGGK_256x64x5x3_nu_scan';
%% Convergence cases kN = (1.6 2.2) nu = (0.01 1.0)
resdir = 'Zpinch_rerun/convcoll_Kn_1.6_200x32x3x2_mu_0.01';

JOBNUMMIN = 00; JOBNUMMAX = 10;
resdir = ['results/',resdir];
run analysis_gyacomo