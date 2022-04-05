% folder = '/misc/gene_results/shearless_cyclone/miller_output_1.0/';
% folder = '/misc/gene_results/shearless_cyclone/miller_output_0.8/';
% folder = '/misc/gene_results/shearless_cyclone/s_alpha_output_1.0/';
folder = '/misc/gene_results/shearless_cyclone/s_alpha_output_0.8/';
% folder = '/misc/gene_results/HP_fig_2b_mu_5e-2/';
% folder = '/misc/gene_results/HP_fig_2c_mu_5e-2/';
gene_data = load_gene_data(folder);
gene_data = rotate_c_plane_nxnky_to_nkxny(gene_data);
if 1
%% Space time diagramm (fig 11 Ivanov 2020)
TAVG_0 = 500; TAVG_1 = 700; % Averaging times duration
% chose your field to plot in spacetime diag (uzf,szf,Gx)
field = 'phi';
compz = 'avg';
nmvm  = 1;
fig = plot_radial_transport_and_spacetime(gene_data,TAVG_0,TAVG_1,field,nmvm,compz);
% save_figure(data,fig)
end

if 0
%% 2D snapshots
% Options
options.INTERP    = 0;
options.POLARPLOT = 0;
options.AXISEQUAL = 1;
options.NAME      = '\phi';
% options.NAME      = 'v_y';
% options.NAME      = 'T_i';
% options.NAME      = '\Gamma_x';
% options.NAME      = 'k^2n_e';
options.PLAN      = 'xy';
% options.NAME      ='f_e';
% options.PLAN      = 'sx';
options.COMP      = 9;
options.TIME      = [100 200 400];
data.a = data.EPS * 2000;
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
options.PLAN      = 'xy';
% options.NAME      = 'f_e';
% options.PLAN      = 'sx';
options.COMP      = 9;
options.TIME      = 0:700;
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
