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
outfile = 'shearless_cyclone/64x32x160x5x3_CBC_Npol_10_kine';
%% CBC
% outfile ='shearless_cyclone/64x32x16x5x3_CBC_080';
% outfile ='shearless_cyclone/64x32x16x5x3_CBC_scan/64x32x16x5x3_CBC_100';
% outfile ='shearless_cyclone/64x32x16x5x3_CBC_120';

% outfile ='shearless_cyclone/64x32x16x9x5_CBC_080';
% outfile ='shearless_cyclone/64x32x16x9x5_CBC_100';
% outfile ='shearless_cyclone/64x32x16x9x5_CBC_120';

% outfile = 'shearless_cyclone/64x32x16x5x3_CBC_CO/64x32x16x5x3_CBC_LRGK';

run analysis_HeLaZ
