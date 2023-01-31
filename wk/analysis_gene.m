gyacomodir = '/home/ahoffman/gyacomo/';
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add

% folder = '/misc/gene_results/shearless_cyclone/miller_output_1.0/';
% folder = '/misc/gene_results/shearless_cyclone/miller_output_0.8/';
% folder = '/misc/gene_results/shearless_cyclone/s_alpha_output_1.0/';
% folder = '/misc/gene_results/shearless_cyclone/rm_corrections_HF/';
% folder = '/misc/gene_results/shearless_cyclone/linear_s_alpha_CBC_100/';
% folder = '/misc/gene_results/shearless_cyclone/s_alpha_output_0.5/';
% folder = '/misc/gene_results/shearless_cyclone/LD_s_alpha_output_1.0/';
% folder = '/misc/gene_results/shearless_cyclone/LD_s_alpha_output_0.8/';
% folder = '/misc/gene_results/Z-pinch/HP_fig_2a_mu_1e-2/';
% folder = '/misc/gene_results/Z-pinch/HP_fig_2a_gyroLES/';
% folder = '/misc/gene_results/Z-pinch/HP_fig_2b_mu_5e-2/';
% folder = '/misc/gene_results/Z-pinch/HP_fig_2c_mu_5e-2/';
% folder = '/misc/gene_results/Z-pinch/HP_fig_2c_gyroLES/';
% folder = '/misc/gene_results/Z-pinch/HP_fig_2c_mu_1e-2_muv_1e-1/';
% folder = '/misc/gene_results/LD_zpinch_1.6/';
% folder = '/misc/gene_results/ZP_HP_kn_1.6_nuv_3.2/';
% folder = '/misc/gene_results/ZP_HP_kn_1.6_nuv_3.2/';
% folder = '/misc/gene_results/Z-pinch/ZP_HP_kn_1.6_HRES/';
% folder = '/misc/gene_results/ZP_kn_2.5_large_box/';
% folder = '/misc/gene_results/Z-pinch/kN_2.0_HD_transport_spectrum_01/';
% folder = '/misc/gene_results/Z-pinch/kN_2.5_HD_transport_spectrum_01/';

% folder = '/misc/gene_results/CBC/128x64x16x24x12/';
% folder = '/misc/gene_results/CBC/196x96x20x32x16_01/';
% folder = '/misc/gene_results/CBC/256x96x24x32x12/';
% folder = '/misc/gene_results/CBC/128x64x16x6x4/';
% folder = '/misc/gene_results/CBC/KT_5.3_128x64x16x24x12_01/';
% folder = '/misc/gene_results/CBC/KT_4.5_128x64x16x24x12_01/';
% folder = '/misc/gene_results/CBC/KT_9_128x64x16x24x12/';
% folder = '/misc/gene_results/CBC/KT_13_large_box_128x64x16x24x12/';
% folder = '/misc/gene_results/CBC/Lapillone_Fig6/';
% folder = '/misc/gene_results/Z-pinch/HP_kN_1.6_adapt_mu_01/';
% folder = '/misc/gene_results/miller/';
% folder = '/misc/gene_results/CBC/KT_4.5_192x96x24x30x16_00/';
% folder = '/misc/gene_results/CBC/KT_5.0_192x96x24x30x16_00/';
% folder = '/misc/gene_results/CBC/KT_5.3_192x96x24x30x16_00/';
% folder = '/misc/gene_results/CBC/196x96x20x32x16_01/';
% folder = '/misc/gene_results/linear_miller_CBC/hdf5_miller_s0_adiabe/';
% folder = '/misc/gene_results/linear_miller_CBC/hdf5_salpha_s0_adiabe/';

%Paper 2
% folder = '/misc/gene_results/CBC/KT_6.96_64x32x32x24x12_Nexc_5/';
% folder = '/misc/gene_results/CBC/KT_6.96_128x64x24x8x4_Nexc_5_00/';
% folder = '/misc/gene_results/CBC/KT_6.96_128x64x24x16x8_Nexc_5_00/';
% folder = '/misc/gene_results/CBC/KT_6.96_128x64x24x32x16_Nexc_5_00/';
% folder = '/misc/gene_results/CBC/KT_6.96_128x64x24x32x16_Nexc_5_01/';

% folder = '/misc/gene_results/CBC/KT_5.3_128x64x24x32x16_Nexc_5_00/';
% folder = '/misc/gene_results/CBC/KT_5.3_128x64x24x32x16_Nexc_5_01/';
% folder = '/misc/gene_results/CBC/new_sim/KT_5.3_128x64x24x16x8_Nexc_5/';
% folder = '/misc/gene_results/CBC/new_sim/KT_5.3_128x64x24x8x4_Nexc_5/';
folder = '/misc/gene_results/CBC/new_sim/KT_6.96_128x64x24x8x4_Nexc_5_smallvbox/';
% folder = '/misc/gene_results/CBC/new_sim/KT_6.96_128x64x24x16x8_Nexc_5_largexbox/';

gene_data = load_gene_data(folder);
gene_data.FIGDIR = folder;
gene_data = invert_kxky_to_kykx_gene_results(gene_data);
gene_data.Pmaxi = gene_data.Nvp-1;
gene_data.Jmaxi = gene_data.Nmu-1;
gene_data.CODENAME = 'GENE';

%% Dashboard (Compilation of main plots of the sim)
dashboard(gene_data);

%% Separated plot routines
if 0
%% Space time diagramm (fig 11 Ivanov 2020)
options.TAVG_0   = 0.1*gene_data.Ts3D(end);
options.TAVG_1   = gene_data.Ts3D(end); % Averaging times duration
options.NMVA     = 1;              % Moving average for time traces
options.ST_FIELD = '\phi';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x, Q_x)
options.INTERP   = 1;
options.NCUT     = 4;              % Number of cuts for averaging and error estimation
options.RESOLUTION = 256;
fig = plot_radial_transport_and_spacetime(gene_data,options,'GENE');
save_figure(gene_data,fig,'.png')
end

if 0
%% statistical transport averaging
options.T = [100 500];
fig = statistical_transport_averaging(gene_data,options);
end

if 0
%% fields snapshots
% Options
options.INTERP    = 0;
options.POLARPLOT = 0;
options.AXISEQUAL = 0;
options.NORMALIZE = 0;
% options.NAME      = '\phi';
% options.NAME      = '\psi';
% options.NAME      = '\omega_z';
options.NAME      = 'n_i';
% options.NAME      = 'n_i-n_e';
% options.NAME      = '\phi^{NZ}';
% options.NAME      = 'N_i^{00}';
% options.NAME      = 'N_i^{00}-N_e^{00}';
% options.NAME      = 's_{Ex}';
% options.NAME      = 'Q_x';
% options.NAME      = 'k^2n_e';
options.PLAN      = 'xy';
options.COMP      = 1;
options.TIME      = [0];
options.RESOLUTION = 256;

data.a = data.EPS * 2e3;
fig = photomaton(gene_data,options);
% save_figure(data,fig)
end

if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
options.INTERP    = 1;
options.POLARPLOT = 0;
options.NAME      = '\phi';
% options.NAME      = 'v_{Ey}';
% options.NAME      = 'G_x';
% options.NAME      = 'n_i';
options.PLAN      = 'xy';
% options.NAME      = 'f_e';
% options.PLAN      = 'sx';
options.COMP      = 'avg';
options.TIME      = gene_data.Ts3D;
options.RESOLUTION= 256;
gene_data.a = data.EPS * 2000;
create_film(gene_data,options,'.gif')
end

if 0
%% Geometry
names = {'$g^{xx}$','$g^{xy}$','$g^{xz}$','$g^{yy}$','$g^{yz}$','$g^{zz}$',...
         '$B_0$','$\partial_x B_0$','$\partial_y B_0$','$\partial_z B_0$',...
         '$J$','$R$','$\phi$','$Z$','$\partial_R x$','$\partial_Z x$'};
figure;
subplot(311)
    for i = 1:6
    plot(gene_data.z, gene_data.geo_arrays(:,i),'DisplayName',names{i}); hold on;
    end
    xlim([min(gene_data.z),max(gene_data.z)]); legend('show'); title('GENE geometry');

subplot(312)
    for i = 7:10
    plot(gene_data.z, gene_data.geo_arrays(:,i),'DisplayName',names{i}); hold on;
    end
    xlim([min(gene_data.z),max(gene_data.z)]); legend('show');

subplot(313)
    for i = [11 12 14]
    plot(gene_data.z, gene_data.geo_arrays(:,i),'DisplayName',names{i}); hold on;
    end
    xlim([min(gene_data.z),max(gene_data.z)]); legend('show');

end

if 0
%% Show f_i(vpar,mu)
options.T         = [0.5 1]*gene_data.Ts3D(end);
options.SPECIES   = 'i';
% options.PLT_FCT = 'contour';
% options.PLT_FCT = 'contourf';
% options.PLT_FCT = 'surf';
options.PLT_FCT = 'surfvv';
options.non_adiab = 0;
options.RMS       = 1; % Root mean square i.e. sqrt(sum_k|f_k|^2) as in Gene
options.folder  = gene_data.folder;
options.iz      = 'avg';
options.FIELD   = '<f_>';
options.SPAR    = linspace(-3,3,32);
options.XPERP   = linspace( 0,sqrt(6),16).^2;
options.ONED    = 1;
plot_fa_gene(options);
end

if 0
%% Time averaged spectrum
options.TIME   = [1];
options.NORM   =1;
options.NAME   = '\phi';
% options.NAME      = 'n_i';
% options.NAME   ='\Gamma_x';
options.PLAN   = 'kxky';
options.COMPZ  = 'avg';
options.OK     = 0;
options.COMPXY = 'avg'; % avg/sum/max/zero/ 2D plot otherwise
options.COMPT  = 'avg';
options.PLOT   = 'semilogy';
fig = spectrum_1D(gene_data,options);
% save_figure(data,fig)
end

if 0
%% Mode evolution
options.NORMALIZED = 0;
options.K2PLOT = 1;
options.TIME   = 1:700;
options.KX_TW  = [25 55]; %kx Growth rate time window
options.KY_TW  = [0 20];  %ky Growth rate time window
options.NMA    = 1;
options.NMODES = 15;
options.iz     = 'avg';
options.ik     = 1; % sum, max or index
options.fftz.flag = 0;
fig = mode_growth_meter(gene_data,options);
save_figure(gene_data,fig)
end
