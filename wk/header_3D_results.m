% Directory of the code "mypathtoHeLaZ/HeLaZ/"
gyacomodir = '/home/ahoffman/gyacomo/';
% Directory of the simulation (from results)
% if 1% Local results
% resdir ='volcokas/64x32x16x5x3_kin_e_npol_1';

%% Dimits
% resdir ='shearless_cyclone/128x64x16x5x3_Dim_90';
% resdir ='shearless_cyclone/128x64x16x9x5_Dim_scan/128x64x16x9x5_Dim_60';
% resdir ='shearless_cyclone/128x64x16x5x3_Dim_scan/128x64x16x5x3_Dim_70';
% resdir ='shearless_cyclone/64x32x16x5x3_Dim_scan/64x32x16x5x3_Dim_70';
%% AVS
% resdir = 'volcokas/64x32x16x5x3_kin_e_npol_1';
% resdir = 'volcokas/64x32x16x5x3_kin_e_npol_1';
% resdir = 'shearless_cyclone/64x32x80x5x3_CBC_Npol_5_kine';
% resdir = 'shearless_cyclone/96x32x160x5x3_CBC_Npol_10_kine';
% resdir = 'shearless_cyclone/64x32x160x5x3_CBC_Npol_10_kine';
% resdir = 'shearless_cyclone/96x32x160x5x3_CBC_Npol_10_kine';
%% shearless CBC
% resdir ='shearless_cyclone/64x32x16x5x3_CBC_080';
% resdir ='shearless_cyclone/64x32x16x5x3_CBC_scan/64x32x16x5x3_CBC_100';
% resdir ='shearless_cyclone/64x32x16x5x3_CBC_120';

% resdir ='shearless_cyclone/64x32x16x9x5_CBC_080';
% resdir ='shearless_cyclone/64x32x16x9x5_CBC_100';
% resdir ='shearless_cyclone/64x32x16x9x5_CBC_120';

% resdir = 'shearless_cyclone/64x32x16x5x3_CBC_CO/64x32x16x5x3_CBC_LRGK';


%% CBC
% resdir = 'CBC/64x32x16x5x3';
% resdir = 'CBC/64x128x16x5x3';
% resdir = 'CBC/128x64x16x5x3';
% resdir = 'CBC/96x96x16x3x2_Nexc_6';
% resdir = 'CBC/128x96x16x3x2';
% resdir = 'CBC/192x96x16x3x2';
resdir = 'CBC/192x96x24x13x7';
% resdir = 'CBC/kT_11_128x64x16x5x3';
% resdir = 'CBC/kT_9_256x128x16x3x2';
% resdir = 'CBC/kT_4.5_128x64x16x13x3';
% resdir = 'CBC/kT_4.5_192x96x24x13x7';
% resdir = 'CBC/kT_4.5_128x64x16x13x7';
% resdir = 'CBC/kT_4.5_128x96x24x15x5';
% resdir = 'CBC/kT_5.3_192x96x24x13x7';
% resdir = 'CBC/kT_13_large_box_128x64x16x5x3';
% resdir = 'CBC/kT_11_96x64x16x5x3_ky_0.02';

% resdir = 'CBC/kT_scan_128x64x16x5x3';
% resdir = 'CBC/kT_scan_192x96x16x3x2';
% resdir = 'CBC/kT_13_96x96x16x3x2_Nexc_6';
% resdir = 'dbg/nexc_dbg';
% resdir = 'CBC/NM_F4_kT_4.5_192x64x24x6x4';

% resdir = 'CBC_Ke_EM/192x96x24x5x3';
% resdir = 'CBC_Ke_EM/96x48x16x5x3';
% resdir = 'CBC_Ke_EM/minimal_res';
%% KBM
% resdir = 'NL_KBM/192x64x24x5x3';
%% Linear CBC
% resdir = 'linear_CBC/20x2x32_21x11_Lx_62.8319_Ly_31.4159_q0_1.4_e_0.18_s_0.8_kN_2.22_kT_5.3_nu_1e-02_DGDK_adiabe';
% resdir = 'testcases/miller_example';
% resdir = 'Miller/128x256x3x2_CBC_dt_5e-3';
resdir = ['results/',resdir];
JOBNUMMIN = 00; JOBNUMMAX = 10;
run analysis_gyacomo
