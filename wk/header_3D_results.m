% Directory of the code "mypathtoHeLaZ/HeLaZ/"
helazdir = '/home/ahoffman/HeLaZ/';
% Directory of the simulation (from results)
% if 1% Local results
outfile ='';

%% Dimits
% outfile ='shearless_cyclone/128x64x16x5x3_Dim_90';
outfile ='shearless_cyclone/64x32x16x5x3_Dim_50';
%% AVS
% outfile = 'volcokas/64x32x16x5x3_kin_e_npol_1';

%% Bechmark
% outfile ='shearless_cyclone/64x32x16x5x3_CBC_080';
% outfile ='shearless_cyclone/64x32x16x5x3_CBC_100';
% outfile ='shearless_cyclone/64x32x16x5x3_CBC_120';

% outfile ='shearless_cyclone/64x32x16x9x5_CBC_080';
% outfile ='shearless_cyclone/64x32x16x9x5_CBC_100';
% outfile ='shearless_cyclone/64x32x16x9x5_CBC_120';

run analysis_HeLaZ
