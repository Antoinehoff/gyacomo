%% Directory of the simulation
% if 1% Local results
outfile ='';
outfile ='';
outfile ='';
outfile ='';
% outfile ='debug/ppj_init';
%% nu = 5e-1
% Sugama
% outfile ='Hallenbert_nu_5e-01/200x32_5x3_L_120_kN_1.5_kT_0.375_nu_5e-01_SGGK';% also in 7x4
% outfile ='Hallenbert_nu_5e-01/200x32_5x3_L_120_kN_1.6_kT_0.4_nu_5e-01_SGGK';
% outfile ='Hallenbert_nu_5e-01/200x32_7x4_L_120_kN_1.7_kT_0.425_nu_5e-01_SGGK';% also in 7x4
% outfile ='Hallenbert_nu_5e-01/200x32_5x3_L_120_kN_1.8_kT_0.45_nu_5e-01_SGGK';
% outfile ='Hallenbert_nu_5e-01/200x32_5x3_L_120_kN_1.9_kT_0.475_nu_5e-01_SGGK';%also in 7x4

%% nu = 1e-1
% Landau
% outfile ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.5_kT_0.375_nu_1e-01_LDGK';
% outfile ='Hallenbert_nu_1e-01/150x50_5x3_L_120_kN_1.6_kT_0.4_nu_1e-01_LDGK';
% outfile ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.6_kT_0.4_nu_1e-01_LDGK';
% outfile ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.7_kT_0.425_nu_1e-01_LDGK';
% outfile ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.8_kT_0.45_nu_1e-01_LDGK';
% outfile ='Hallenbert_nu_1e-01/150x50_5x3_L_120_kN_1.8_kT_0.45_nu_1e-01_LDGK';
% outfile ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.9_kT_0.475_nu_1e-01_LDGK';

% Sugama
% outfile ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.5_kT_0.375_nu_1e-01_SGGK';
% outfile ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.7_kT_0.425_nu_1e-01_SGGK';
% outfile ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.8_kT_0.45_nu_1e-01_SGGK';
% outfile ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.9_kT_0.475_nu_1e-01_SGGK';

% Dougherty
% outfile ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.5_kT_0.375_nu_1e-01_DGGK';
% outfile ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.6_kT_0.4_nu_1e-01_DGGK';
% outfile ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.7_kT_0.425_nu_1e-01_DGGK';
% outfile ='Hallenbert_nu_1e-01/200x32_5x3_L_120_kN_1.8_kT_0.45_nu_1e-01_DGGK';

%% nu = 5e-2
% outfile ='Hallenbert_nu_5e-02/200x32_11x6_L_120_kN_1.8_kT_0.45_nu_5e-02_SGGK';%For GENE benchmark % to analyse (added HD)
% outfile ='Hallenbert_nu_5e-02/200x32_11x6_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_SGGK';

% testing various NL closures
% outfile ='Hallenbert_nu_5e-02/200x32_7x4_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_SGGK';

% outfile ='Hallenbert_nu_5e-02/200x32_5x3_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_SGDK';
% outfile ='Hallenbert_nu_5e-02/200x32_11x6_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_SGDK';

% outfile ='Hallenbert_nu_5e-02/256x64_5x3_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_SGDK';
% outfile ='Hallenbert_nu_5e-02/256x64_11x6_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_SGDK';
% outfile ='Hallenbert_nu_5e-02/200x32_21x3_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_SGDK';
% outfile ='Hallenbert_nu_5e-02/200x32_17x9_L_120_kN_1.8_kT_0.45_nu_5e-02_SGDK';

% outfile ='Hallenbert_nu_5e-02/200x32_5x3_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_DGGK';
% outfile ='Hallenbert_nu_5e-02/200x32_11x6_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_DGGK';
% outfile ='Hallenbert_nu_5e-02/128x32_5x3_Lx_120_Ly_60_kN_1.8_kT_0.45_nu_5e-02_FCGK';

%% nu = 1e-2
% Landau
% outfile ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.5_kT_0.375_nu_1e-02_LDGK';
% outfile ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.6_kT_0.4_nu_1e-02_LDGK';
% outfile ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.7_kT_0.425_nu_1e-02_LDGK';
% outfile ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.8_kT_0.45_nu_1e-02_LDGK';
% outfile ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.9_kT_0.475_nu_1e-02_LDGK';

% Sugama
% outfile ='kobayashi_2015_fig1/150x150_5x3_L_100_kN_1.4_nu_5e-03_SGGK';
% outfile ='Hallenbert_nu_1e-02/200x32_7x4_L_120_kN_1.5_kT_0.375_nu_1e-02_SGGK';
% outfile ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.5_kT_0.375_nu_1e-02_SGGK';
% outfile ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.6_kT_0.4_nu_1e-02_SGGK';
% outfile ='Hallenbert_nu_1e-02/200x32_7x4_L_120_kN_1.6_kT_0.4_nu_1e-02_SGGK';
% outfile ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.7_kT_0.425_nu_1e-02_SGGK';
% outfile ='Hallenbert_nu_1e-02/300x64_5x3_L_120_kN_1.7_kT_0.425_nu_1e-02_SGGK';
% outfile ='Hallenbert_nu_1e-02/200x32_11x6_L_120_kN_1.7_kT_0.425_nu_1e-02_SGGK';
% outfile ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.8_kT_0.45_nu_1e-02_SGGK'; % To analyse (added HD)
% outfile ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.9_kT_0.475_nu_1e-02_SGGK';
% outfile ='Hallenbert_nu_1e-02/200x32_11x6_L_120_kN_1.9_kT_0.475_nu_1e-02_SGGK';
% Dougherty
% outfile ='Hallenbert_nu_1e-02/200x32_7x4_L_120_kN_1.5_kT_0.375_nu_1e-02_DGGK';
% outfile ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.6_kT_0.4_nu_1e-02_DGGK';
% outfile ='Hallenbert_nu_1e-02/200x32_8x5_L_120_kN_1.6_kT_0.4_nu_0e+00_DGGK';
% outfile ='Hallenbert_nu_1e-02/200x32_5x3_L_120_kN_1.8_kT_0.45_nu_1e-02_DGGK';

%% nu = 5e-3
% outfile ='Hallenbert_nu_5e-03/200x32_5x3_Lx_120_Ly_60_kN_1.8_eta_0.25_nuSG_5e-03_muxy_5e-2';
% outfile ='Hallenbert_nu_5e-03/200x32_5x3_Lx_120_Ly_60_kN_1.8_eta_0.25_nuSG_5e-03_mux_5e-2_muy_6e-1';
% outfile ='Hallenbert_nu_5e-03/200x32_11x6_Lx_120_Ly_60_kN_1.8_eta_0.25_nuSG_5e-03_mux_5e-2_muy_6e-1';
% outfile ='Hallenbert_nu_5e-03/200x32_11x6_Lx_120_Ly_60_kN_1.8_eta_0.25_nuSG_5e-03_muxy_5e-2';

%% nu = 0

% outfile ='Hallenbert_fig2a/200x32_21x11_Lx_120_Ly_60_kN_1.6_eta_0.4_nu_0_muxy_1e-2';
% outfile ='Hallenbert_fig2a/200x32_11x6_Lx_120_Ly_60_kN_1.6_eta_0.4_nu_0_muxy_1e-2';
% outfile ='Hallenbert_fig2a/200x32_5x3_Lx_120_Ly_60_kN_1.6_eta_0.4_nu_0_muxy_1e-2';
% outfile ='Hallenbert_fig2a/200x32_5x3_Lx_120_Ly_60_kN_1.6_eta_0.4_nuDGGK_0.1_muxy_1e-2';
% outfile ='Hallenbert_fig2a/200x32_11x6_Lx_120_Ly_60_kN_1.6_eta_0.4_nuDGGK_0.1_muxy_1e-2';
% outfile ='Hallenbert_fig2a/200x32_5x3_Lx_120_Ly_60_kN_1.6_eta_0.4_nuSGGK_0.1_muxy_1e-2';
% outfile ='Hallenbert_fig2a/200x32_11x6_Lx_120_Ly_60_kN_1.6_eta_0.4_nuSGGK_0.1_muxy_1e-2';
% outfile ='Hallenbert_fig2a/200x32_5x3_Lx_120_Ly_60_kN_1.6_eta_0.4_nuLDGK_0.1_muxy_1e-2';
% outfile ='Hallenbert_fig2a/200x32_5x3_Lx_120_Ly_60_kN_1.6_eta_0.4_nuLRGK_0.1_muxy_1e-2';


% outfile ='Hallenbert_fig2b/200x32_11x6_Lx_240_Ly_120_kN_2.5_eta_0.25_nu_0_muxy_1e-1';
% outfile ='Hallenbert_fig2b/200x32_5x3_Lx_240_Ly_120_kN_2.5_eta_0.25_nu_0_muxy_1e-1';
% outfile ='Hallenbert_fig2b/200x32_5x3_Lx_240_Ly_120_kN_2.5_eta_0.25_nuSGGK_0.1_muxy_1e-1';
% outfile ='Hallenbert_fig2b/200x32_5x3_Lx_240_Ly_120_kN_2.5_eta_0.25_nuDGGK_0.1_muxy_1e-1';
% outfile ='Hallenbert_fig2b/200x32_5x3_Lx_240_Ly_120_kN_2.5_eta_0.25_nuLDGK_0.1_muxy_1e-1';
% outfile ='Hallenbert_fig2b/200x32_5x3_Lx_240_Ly_120_kN_2.5_eta_0.25_nuLRGK_0.1_muxy_1e-1';

% outfile ='Hallenbert_fig2c/200x32_11x6_Lx_120_Ly_60_kN_2.0_eta_0.25_nu_0_muxy_5e-2';
% outfile ='Hallenbert_fig2c/200x32_5x3_Lx_120_Ly_60_kN_2.0_eta_0.25_nu_0_muxy_5e-2';
% outfile ='Hallenbert_fig2c/200x32_5x3_Lx_120_Ly_60_kN_2.0_eta_0.25_nuSGGK_0.1_muxy_5e-2';
% outfile ='Hallenbert_fig2c/200x32_11x6_Lx_120_Ly_60_kN_2.0_eta_0.25_nuSGGK_0.1_muxy_5e-2';
% outfile ='Hallenbert_fig2c/200x32_5x3_Lx_120_Ly_60_kN_2.0_eta_0.25_nuDGGK_0.1_muxy_5e-2';
% outfile ='Hallenbert_fig2c/200x32_11x6_Lx_120_Ly_60_kN_2.0_eta_0.25_nuDGGK_0.1_muxy_5e-2';
% outfile ='Hallenbert_fig2c/200x32_5x3_Lx_120_Ly_60_kN_2.0_eta_0.25_nuLDGK_0.1_muxy_5e-2';
% outfile ='Hallenbert_fig2c/200x32_5x3_Lx_120_Ly_60_kN_2.0_eta_0.25_nuLRGK_0.1_muxy_5e-2';
% outfile ='Hallenbert_fig2c/200x32_9x5_Lx_120_Ly_60_kN_2.0_eta_0.25_nuLRGK_0.1_muxy_5e-2';

%% Transport scan
% outfile = 'nu_0.1_transport_scan/colless_kn_1.7_to_2.0';
% outfile = 'nu_0.1_transport_scan/colless_kn_2.1_to_2.5';


% outfile = 'nu_0.1_transport_scan/DG_kn_1.8_to_2.1';
% outfile = 'nu_0.1_transport_scan/DG_kn_2.2_to_2.5';

% outfile = 'nu_0.1_transport_scan/SG_kn_1.7_to_2.0';
% outfile = 'nu_0.1_transport_scan/SG_10x5_conv_test';
% outfile = 'nu_0.1_transport_scan/SG_kn_2.2_to_2.5';

% outfile = 'nu_0.1_transport_scan/LD_kn_2.0_to_2.5';

% outfile = 'nu_0.1_transport_scan/LR_kn_1.7_to_2.0';
% outfile = 'nu_0.1_transport_scan/LR_kn_2.1_to_2.5';

% outfile = 'nu_0.1_transport_scan/colless_kn_2.2_Lx1.5';
% outfile = 'nu_0.1_transport_scan/colless_kn_2.2_HD';

% outfile = 'nu_0.1_transport_scan/colless_kn_1.6_HD';

% outfile = 'nu_0.1_transport_scan/large_box_kN_2.1_nu_0.1';

% outfile = 'predator_prey/SG_Kn_1.7_nu_0.01';

outfile = 'ZF_damping_DK/SG_4x2_nu_0.1';
% else% Marconi results
% outfile ='';
% outfile ='';
% outfile ='';
% outfile ='';
% outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/simulation_A_new/300x300_5x3_L_120_kN_1.6667_nu_1e-01_SGGK/out.txt';
% % outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/simulation_A/300x150_L_120_P_8_J_4_eta_0.6_nu_1e-01_SGGK_mu_0e+00/out.txt';
% % BASIC.RESDIR      = ['../',outfile(46:end-8),'/'];
% MISCDIR = ['/misc/HeLaZ_outputs/',outfile(46:end-8),'/'];
% end

analysis_3D