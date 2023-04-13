% Directory of the code "mypathtogyacomo/gyacomo/"
% Partition of the computer where the data have to be searched
% PARTITION  = '/misc/gyacomo23_outputs/';
PARTITION = '/home/ahoffman/gyacomo/';
%% CBC 
% resdir = 'paper_2_GYAC23/CBC/5x3x128x64x24_dp';
% resdir = 'paper_2_GYAC23/CBC/7x4x128x64x24_dp';
% resdir = 'paper_2_GYAC23/CBC/9x5x128x64x24_dp';
% resdir = 'paper_2_GYAC23/CBC/9x5x192x96x32_dp';
% resdir = 'paper_2_GYAC23/CBC/11x6x128x64x24_dp';
% resdir = 'paper_2_GYAC23/CBC/11x6x128x64x24_dp';
% resdir = 'paper_2_GYAC23/CBC/21x11x128x64x24_dp';

%% tests single vs double precision
% resdir = 'paper_2_GYAC23/precision_study/5x3x128x64x24';
% resdir = 'paper_2_GYAC23/precision_study/5x3x128x64x24_dp';
% resdir = 'paper_2_GYAC23/precision_study/5x3x128x64x24_sp';
% resdir = 'paper_2_GYAC23/precision_study/5x3x128x64x24_sp_clos_1';
% resdir = 'paper_2_GYAC23/precision_study/3x2x128x64x24_sp_muz_2.0';
% resdir = 'paper_2_GYAC23/precision_study/test_3x2x128x64x24_sp_muz_2.0';
% resdir = 'paper_2_GYAC23/precision_study/3x2x128x64x24_sp_clos_1';

%% Marconi results
% resdir = 'paper_2_GYAC23/collisionless/kT_5.3/5x3x128x64x24_dp_muz_2.0_muxy_0';
% resdir = 'paper_2_GYAC23/collisionless/kT_5.3/5x3x128x64x24_dp_SG';
% resdir = 'paper_2_GYAC23/collisionless/kT_5.3/5x3x128x64x24_dp_muz_2.0_full_NL';
% resdir = 'paper_2_GYAC23/collisionless/kT_5.3/7x4x128x64x24_dp';
% resdir = 'paper_2_GYAC23/collisionless/kT_5.3/9x5x128x64x24_dp';
% resdir = 'paper_2_GYAC23/collisionless/kT_5.3/11x6x128x64x24_dp';

% resdir = 'paper_2_GYAC23/collisionless/CBC/5x3x128x64x24_dp';
% resdir = 'paper_2_GYAC23/collisionless/CBC/7x4x128x64x24_dp';
% resdir = 'paper_2_GYAC23/collisionless/CBC/9x5x128x64x24_dp';
% resdir = 'paper_2_GYAC23/collisionless/CBC/11x6x128x64x24_dp';
% resdir = 'paper_2_GYAC23/collisionless/CBC/9x5x192x96x32_dp';

% resdir = 'paper_2_GYAC23/collisionless/kT_scan_nu_1e-3/5x3x128x64x24_dp';
% resdir = 'paper_2_GYAC23/collisionless/kT_scan_nu_1e-3/7x4x128x64x24_dp';
% resdir = 'paper_2_GYAC23/collisionless/kT_scan_nu_1e-3/9x5x128x64x24_dp';

% resdir = 'paper_2_GYAC23/collision_study/nuDGGK_scan_kT_5.3/5x3x128x64x24_dp';
% resdir = 'paper_2_GYAC23/collision_study/nuDGGK_scan_kT_5.3/9x5x128x64x24_dp';

%% low precision 3D ITG
% resdir = 'results/paper_2_GYAC23/3x2x64x48x16/CBC_3x2x64x48x16_CLOS_1';
% resdir = 'results/paper_2_GYAC23/3x2x64x48x16/kT_0.0';
% resdir = 'results/paper_2_GYAC23/3x2x64x48x16/kT_3.0';
% resdir = 'results/paper_2_GYAC23/3x2x64x48x16/kT_3.5';
% resdir = 'results/paper_2_GYAC23/3x2x64x48x16/kT_4.0';
% resdir = 'results/paper_2_GYAC23/3x2x64x48x16/kT_4.5';
% resdir = 'results/paper_2_GYAC23/3x2x64x48x16/kT_5.3';
% resdir = 'results/paper_2_GYAC23/3x2x64x48x16/CBC';

% resdir = 'results/paper_2_GYAC23/5x2x64x48x16/kT_3.5';
% resdir = 'results/paper_2_GYAC23/5x2x64x48x16/kT_4.0';
% resdir = 'results/paper_2_GYAC23/5x2x64x48x16/kT_4.5';
% resdir = 'results/paper_2_GYAC23/5x2x64x48x16/kT_5.3';
% resdir = 'results/paper_2_GYAC23/5x2x64x48x16/CBC';

% resdir = 'results/paper_2_GYAC23/9x2x64x48x16/kT_3.5';
% resdir = 'results/paper_2_GYAC23/9x2x64x48x16/kT_4.0';
% resdir = 'results/paper_2_GYAC23/9x2x64x48x16/kT_4.5';
% resdir = 'results/paper_2_GYAC23/9x2x64x48x16/kT_5.3';
% resdir = 'results/paper_2_GYAC23/9x2x64x48x16/CBC';

% resdir = 'results/paper_2_GYAC23/11x2x64x48x16/kT_3.5';
% resdir = 'results/paper_2_GYAC23/11x2x64x48x16/kT_4.0';
% resdir = 'results/paper_2_GYAC23/11x2x64x48x16/kT_4.5';
% resdir = 'results/paper_2_GYAC23/11x2x64x48x16/kT_5.3';
% resdir = 'results/paper_2_GYAC23/11x2x64x48x16/CBC';

%% Box size effect on CBC
% resdir = 'results/paper_2_GYAC23/7x4x128x64x24/CBC_L120';
resdir = 'results/paper_2_GYAC23/7x4x128x64x24/CBC_L180';

%% testcases
% resdir = 'testcases/ITG_zpinch';
% resdir = 'testcases/zpinch_example/results_trunc';
% resdir = 'testcases/zpinch_example/results_maxd=2';
% resdir = 'testcases/DLRA_zpinch/base_case';
% resdir = 'testcases/DLRA_zpinch/nsv_filter_2';
% resdir = 'testcases/DLRA_zpinch/nsv_filter_6';
% resdir = 'testcases/cyclone_example';

 %%
J0 = 00; J1 = 10;

% Load basic info (grids and time traces)
DATADIR = [PARTITION,resdir,'/'];
data    = {};
data    = compile_results_low_mem(data,DATADIR,J0,J1);

if 1
%% Plot transport and phi radial profile
[data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');

options.TAVG_0   = 100;
options.TAVG_1   = 1000;
options.NCUT     = 5;              % Number of cuts for averaging and error estimation
options.NMVA     = 1;              % Moving average for time traces
% options.ST_FIELD = '\Gamma_x';   % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.ST_FIELD = '\phi';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.INTERP   = 0;
options.RESOLUTION = 256;
fig = plot_radial_transport_and_spacetime(data,options);
end

if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
options.INTERP    = 1;
options.POLARPLOT = 0;
options.NAME      = '\phi';
% options.NAME      = '\omega_z';
% options.NAME     = 'N_i^{00}';
% options.NAME      = 's_{Ey}';
% options.NAME      = 'n_i^{NZ}';
% options.NAME      = 'Q_x';
% options.NAME      = 'n_i';
% options.NAME      = 'n_i-n_e';
options.PLAN      = 'xy';
% options.NAME      = 'f_i';
% options.PLAN      = 'sx';
options.COMP      = 'avg';
% options.TIME      = data.Ts5D(end-30:end);
options.TIME      =  data.Ts3D;
% options.TIME      = [0:1500];
data.EPS          = 0.1;
data.a = data.EPS * 2000;
options.RESOLUTION = 256;
create_film(data,options,'.gif')
end

if 0
%% fields snapshots
% Options
[data.Na00, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'Na00');
data.Ni00 = reshape(data.Na00(1,:,:,:,:),data.grids.Nky,data.grids.Nkx,data.grids.Nz,numel(data.Ts3D));

options.INTERP    = 0;
options.POLARPLOT = 0;
options.AXISEQUAL = 0;
options.NORMALIZE = 0;
options.NAME      = 'N_i^{00}';
% options.NAME      = '\phi';
options.PLAN      = 'xy';
options.COMP      = 'avg';
options.TIME      = [800 900 1000];
options.RESOLUTION = 256;
fig = photomaton(data,options);
% save_figure(data,fig)
end
if 0
%% Performance profiler
profiler(data)
end

if 0
%% Hermite-Laguerre spectrum
% [data.Napjz, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'Napjz');
[data.Napjz, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'Nipjz');
options.ST         = 0;
options.NORMALIZED = 0;
options.LOGSCALE   = 1;
options.FILTER     = 0; %filter the 50% time-average of the spectrum from
fig = show_moments_spectrum(data,options);
end

if 0
%% Mode evolution
[data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');

options.NORMALIZED = 0;
options.TIME   = [000:9000];
options.KX_TW  = [1 80]; %kx Growth rate time window
options.KY_TW  = [0 80];  %ky Growth rate time window
options.NMA    = 1;
options.NMODES = 800;
options.iz     = 'avg'; % avg or index
options.ik     = 1; % sum, max or index
options.fftz.flag = 0;
fig = mode_growth_meter(data,options);
% save_figure(data,fig,'.png')
end


if 0
%% Study singular values
[data.SV_ky_pj, data.Ts2D] = compile_results_2D(DATADIR,J0,J1,'sv_ky_pj');
nSV = data.grids.Np * data.grids.Nj;
colors_ = jet(nSV);
figure
for i = 1:nSV
    sv = squeeze(data.SV_ky_pj(i,:));
    semilogy(data.Ts2D,sv,...
        'color',colors_(i,:),'DisplayName',['SV ',num2str(i)]);hold on
end
legend('show');
end