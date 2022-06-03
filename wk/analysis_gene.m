% folder = '/misc/gene_results/shearless_cyclone/miller_output_1.0/';
% folder = '/misc/gene_results/shearless_cyclone/miller_output_0.8/';
folder = '/misc/gene_results/shearless_cyclone/s_alpha_output_1.0/';
% folder = '/misc/gene_results/shearless_cyclone/linear_s_alpha_CBC_100/';
% folder = '/misc/gene_results/shearless_cyclone/s_alpha_output_0.5/';
% folder = '/misc/gene_results/shearless_cyclone/LD_s_alpha_output_1.0/';
% folder = '/misc/gene_results/shearless_cyclone/LD_s_alpha_output_0.8/';
% folder = '/misc/gene_results/HP_fig_2b_mu_5e-2/';
% folder = '/misc/gene_results/HP_fig_2c_mu_5e-2/';
gene_data = load_gene_data(folder);
gene_data = invert_kxky_to_kykx_gene_results(gene_data);
if 1
%% Space time diagramm (fig 11 Ivanov 2020)
options.TAVG_0   = 0.8*gene_data.Ts3D(end);
options.TAVG_1   = gene_data.Ts3D(end); % Averaging times duration
options.NMVA     = 1;              % Moving average for time traces
options.ST_FIELD = '\phi';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x, Q_x)
options.INTERP   = 1;
gene_data.FIGDIR = folder;;
fig = plot_radial_transport_and_spacetime(gene_data,options);
save_figure(gene_data,fig)
end

if 0
%% 2D snapshots
% Options
options.INTERP    = 0;
options.POLARPLOT = 0;
options.AXISEQUAL = 1;
% options.NAME      = 'Q_x';
options.NAME      = '\phi';
% options.NAME      = 'T_i';
% options.NAME      = '\Gamma_x';
% options.NAME      = 'k^2n_e';
options.PLAN      = 'kxky';
% options.NAME      ='f_e';
% options.PLAN      = 'sx';
options.COMP      = 9;
options.TIME      = [0 50 100 200 300];
gene_data.a = data.EPS * 2000;
fig = photomaton(gene_data,options);
save_figure(gene_data,fig)
end

if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
options.INTERP    = 1;
options.POLARPLOT = 0;
options.NAME      = '\phi';
% options.NAME      = 'v_y';
% options.NAME      = 'n_i^{NZ}';
% options.NAME      = '\Gamma_x';
% options.NAME      = 'n_i';
options.PLAN      = 'xz';
% options.NAME      = 'f_e';
% options.PLAN      = 'sx';
options.COMP      = 'avg';
options.TIME      = 000:300;
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
    for i = 11:16
    plot(gene_data.z, gene_data.geo_arrays(:,i),'DisplayName',names{i}); hold on;
    end
    xlim([min(gene_data.z),max(gene_data.z)]); legend('show');

end

if 0
%% Show f_i(vpar,mu)
options.times   = 200:300;
options.specie  = 'i';
options.PLT_FCT = 'contour';
options.folder  = folder;
options.iz      = 9;
options.FIELD   = '<f_>';
options.ONED    = 0;
% options.FIELD   = 'Q_es';
plot_fa_gene(options);
end

if 0
%% Mode evolution
options.NORMALIZED = 0;
options.K2PLOT = 1;
options.TIME   = 100:200;
options.NMA    = 1;
options.NMODES = 15;
options.iz     = 9;
fig = mode_growth_meter(gene_data,options);
save_figure(gene_data,fig)
end
