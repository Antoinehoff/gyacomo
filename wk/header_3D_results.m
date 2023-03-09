% Directory of the code "mypathtoHeLaZ/HeLaZ/"
gyacomodir = '/home/ahoffman/gyacomo/';
% Partition of the computer where the data have to be searched
PARTITION  = '/misc/gyacomo_outputs/';
% PARTITION  = gyacomodir;

%% CBC results with nuDGDK = 0.05
% low resolution (Dirichlet)
% resdir = 'paper_2_nonlinear/kT_6.96/3x2x128x64x24';
% resdir = 'paper_2_nonlinear/kT_6.96/5x3x128x64x24'; %+ diff study
% resdir = 'paper_2_nonlinear/kT_6.96/5x3x128x64x24_dv4_diff';
% resdir = 'paper_2_nonlinear/kT_6.96/optimal_muz_nu_5x3x128x64x24';
% resdir = 'paper_2_nonlinear/kT_6.96/7x4x128x64x24';
% resdir = 'paper_2_nonlinear/kT_6.96/9x5x128x64x24';
% low resolution (Cyclic)
% resdir = 'paper_2_nonlinear/kT_6.96/3x2x128x64x24_cyclic_BC';
% resdir = 'paper_2_nonlinear/kT_6.96/5x3x128x64x24_cyclic_BC';
% resdir = 'paper_2_nonlinear/kT_6.96/7x4x128x64x24_cyclic_BC';
% resdir = 'paper_2_nonlinear/kT_6.96/9x5x128x64x24_cyclic_BC';
% low resolution (periodic)
% resdir = 'paper_2_nonlinear/kT_6.96/3x2x128x64x24_periodic_BC';
% resdir = 'paper_2_nonlinear/kT_6.96/5x3x128x64x24_periodic_BC';
% resdir = 'paper_2_nonlinear/kT_6.96/7x4x128x64x24_periodic_BC';
% resdir = 'paper_2_nonlinear/kT_6.96/9x5x128x64x24_periodic_BC';
% high resolution (Dirichlet)
% resdir = 'paper_2_nonlinear/kT_6.96/5x3x192x96x24';
% resdir = 'paper_2_nonlinear/kT_6.96/7x4x192x96x24';
% resdir = 'paper_2_nonlinear/kT_6.96/9x5x192x96x24';
% high resolution (Cyclic)
% resdir = 'paper_2_nonlinear/kT_6.96/5x3x192x96x24_cyclic_BC';
% resdir = 'paper_2_nonlinear/kT_6.96/7x4x192x96x24_cyclic_BC';
% resdir = 'paper_2_nonlinear/kT_6.96/9x5x192x96x24_cyclic_BC';
% high resolution (periodic)
% resdir = 'paper_2_nonlinear/kT_6.96/5x3x192x96x24_periodic_BC';
% resdir = 'paper_2_nonlinear/kT_6.96/7x4x192x96x24_periodic_BC';
% resdir = 'paper_2_nonlinear/kT_6.96/9x5x192x96x24_periodic_BC';
% resdir = 'paper_2_nonlinear/kT_6.96/11x6x192x96x24_periodic_BC';

% high z resolution (Dirichlet
% resdir = 'paper_2_nonlinear/kT_6.96/5x3x128x64x24';
% resdir = 'paper_2_nonlinear/kT_6.96/5x3x128x64x32';
% resdir = 'paper_2_nonlinear/kT_6.96/5x3x192x96x64';
% resdir = 'paper_2_nonlinear/kT_6.96/5x3x128x64x96';

%% CBC results with nuDGDK = 0.01
% resdir = 'paper_2_nonlinear/kT_6.96_nu_0.01/5x3x192x96x24';

%% kT=5.3 results
resdir = 'paper_2_nonlinear/kT_5.3/5x3x128x64x64';
% resdir = 'paper_2_nonlinear/kT_5.3/5x3x128x64x24';
% resdir = 'paper_2_nonlinear/kT_5.3/7x4x128x64x24';
% resdir = 'paper_2_nonlinear/kT_5.3/7x4x128x64x24_MUxy_0';
% resdir = 'paper_2_nonlinear/kT_5.3/7x4x128x64x24_NL_-1';
% resdir = 'paper_2_nonlinear/kT_5.3/7x4x128x64x24_nuDG_0.01';
% resdir = 'paper_2_nonlinear/kT_5.3/7x4x128x64x64';
% resdir = 'paper_2_nonlinear/kT_5.3/7x4x192x96x64';
% resdir = 'paper_2_nonlinear/kT_5.3/9x5x128x64x24';
% resdir = 'paper_2_nonlinear/kT_5.3/9x5x128x64x64';
% resdir = 'paper_2_nonlinear/kT_5.3/11x6x128x64x24';
% resdir = 'paper_2_nonlinear/kT_5.3/11x6x128x64x64';

%% Old stuff
% resdir = 'CBC/kT_4.5_128x64x16x13x7/';
% resdir = 'CBC/kT_5.3_192x96x24x13x7/';

%% Miller
% resdir = 'paper_2_nonlinear/kT_4.15_miller/7x4x128x64x24';
% resdir = 'paper_2_nonlinear/kT_4.15_miller/7x4x128x64x32';
% resdir = 'paper_2_nonlinear/kT_4.15_miller/7x4x128x64x64';
% resdir = 'paper_2_nonlinear/kT_4.15_miller/7x5x192x64x24';
% resdir = 'paper_2_nonlinear/kT_4.15_miller/GX_fig4_7x5x192x64x24';

%% Redo poster
% resdir = 'paper_2_nonlinear/kT_6.96_Nexc1/5x3x128x64x24'; %% same as poster
%%
JOBNUMMIN = 00; JOBNUMMAX = 10;
run analysis_gyacomo
