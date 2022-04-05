function [ fig ] = plot_metric( data )

names = {'$g^{xx}$','$g^{xy}$','$g^{xz}$','$g^{yy}$','$g^{yz}$','$g^{zz}$',...
         '$B_0$','$\partial_x B_0$','$\partial_y B_0$','$\partial_z B_0$',...
         '$J$','$R$','$\phi$','$Z$','$\partial_R x$','$\partial_Z x$'};
data.geo_arrays = load([data.localdir,'geometry.dat']);
fig = figure;
subplot(311)
    for i = 1:6
    plot(data.z, data.geo_arrays(:,i),'DisplayName',names{i}); hold on;
    end
    xlim([min(data.z),max(data.z)]); legend('show'); title('MoNoLiT geometry');
    
subplot(312)
    for i = 7:10
    plot(data.z, data.geo_arrays(:,i),'DisplayName',names{i}); hold on;
    end
    xlim([min(data.z),max(data.z)]); legend('show');

subplot(313)
    for i = 11:16
    plot(data.z, data.geo_arrays(:,i),'DisplayName',names{i}); hold on;
    end
    xlim([min(data.z),max(data.z)]); legend('show');
end

