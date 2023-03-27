% Directory of the code "mypathtogyacomo/gyacomo/"
% Partition of the computer where the data have to be searched
PARTITION  = '/misc/gyacomo23_outputs/';
% PARTITION  = gyacomodir;

%% CBC 
resdir = 'paper_2_GYAC23/CBC/7x4x192x96x32';

%%
JOBNUMMIN = 00; JOBNUMMAX = 10;

%%
DATADIR = [PARTITION,resdir,'/'];
data    = {};
data    = compile_results_low_mem(data,DATADIR,JOBNUMMIN,JOBNUMMAX);


%%
figure;
plot(data.Ts0D,data.HFLUX_X);