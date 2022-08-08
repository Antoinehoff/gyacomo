% Directory of the code "mypathtoHeLaZ/HeLaZ/"
helazdir = '/home/ahoffman/HeLaZ/';
% Directory of the simulation (from results)
% if 1% Local results
% outfile ='volcokas/64x32x16x5x3_kin_e_npol_1';

%% Dimits
% outfile ='shearless_cyclone/128x64x16x5x3_Dim_90';
% outfile ='shearless_cyclone/128x64x16x9x5_Dim_scan/128x64x16x9x5_Dim_60';
% outfile ='shearless_cyclone/128x64x16x5x3_Dim_scan/128x64x16x5x3_Dim_70';
% outfile ='shearless_cyclone/64x32x16x5x3_Dim_scan/64x32x16x5x3_Dim_70';
%% AVS
% outfile = 'volcokas/64x32x16x5x3_kin_e_npol_1';
% outfile = 'volcokas/64x32x16x5x3_kin_e_npol_1';
% outfile = 'shearless_cyclone/64x32x80x5x3_CBC_Npol_5_kine';
% outfile = 'shearless_cyclone/96x32x160x5x3_CBC_Npol_10_kine';
% outfile = 'shearless_cyclone/64x32x160x5x3_CBC_Npol_10_kine';
% outfile = 'shearless_cyclone/96x32x160x5x3_CBC_Npol_10_kine';
%% shearless CBC
% outfile ='shearless_cyclone/64x32x16x5x3_CBC_080';
% outfile ='shearless_cyclone/64x32x16x5x3_CBC_scan/64x32x16x5x3_CBC_100';
% outfile ='shearless_cyclone/64x32x16x5x3_CBC_120';

% outfile ='shearless_cyclone/64x32x16x9x5_CBC_080';
% outfile ='shearless_cyclone/64x32x16x9x5_CBC_100';
% outfile ='shearless_cyclone/64x32x16x9x5_CBC_120';

% outfile = 'shearless_cyclone/64x32x16x5x3_CBC_CO/64x32x16x5x3_CBC_LRGK';
%% ZPINCH
% outfile ='Zpinch_rerun/Kn_2.5_200x48x5x3';
% outfile ='Zpinch_rerun/Kn_2.5_256x128x5x3';
% outfile ='Zpinch_rerun/Kn_2.5_312x196x5x3_Lx_400_Ly_200';
% outfile ='Zpinch_rerun/Kn_2.5_256x64x5x3';
% outfile ='Zpinch_rerun/Kn_2.0_200x48x9x5_large_box';
% outfile ='Zpinch_rerun/Kn_2.0_256x64x9x5_Lx_240_Ly_120';
% outfile ='Zpinch_rerun/Kn_1.6_256x128x7x4';
% outfile ='Zpinch_rerun/Kn_1.6_200x48x11x6';
% outfile ='Zpinch_rerun/Kn_1.6_256x128x21x11';

%% CBC
% outfile = 'CBC/64x32x16x5x3';
% outfile = 'CBC/64x128x16x5x3';
% outfile = 'CBC/128x64x16x5x3';
% outfile = 'CBC/128x96x16x3x2_Nexc_6';
% outfile = 'CBC/192x96x16x3x2';
% outfile = 'CBC/192x96x24x13x7';
% outfile = 'CBC/kT_11_128x64x16x5x3';
% outfile = 'CBC/kT_9_256x128x16x3x2';
% outfile = 'CBC/kT_4.5_128x64x16x13x2';
% outfile = 'CBC/kT_4.5_192x96x24x13x7';
% outfile = 'CBC/kT_5.3_192x96x24x13x7';
% outfile = 'CBC/kT_13_large_box_128x64x16x5x3';
% outfile = 'CBC/kT_13_96x96x16x3x2_Nexc_6';
% outfile = 'CBC/kT_11_96x64x16x5x3_ky_0.02';

% outfile = 'CBC/kT_scan_128x64x16x5x3';
% outfile = 'CBC/kT_scan_192x96x16x3x2';

outfile = 'CBC/kT_13_96x96x16x3x2_Nexc_6';
%% Linear CBC
% outfile = 'linear_CBC/20x2x32_21x11_Lx_62.8319_Ly_31.4159_q0_1.4_e_0.18_s_0.8_kN_2.22_kT_5.3_nu_1e-02_DGDK_adiabe';

JOBNUMMIN = 00; JOBNUMMAX = 20;
run analysis_HeLaZ
