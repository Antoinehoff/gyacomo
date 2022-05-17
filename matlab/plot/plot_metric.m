function [ fig ] = plot_metric( data )

names = {'Jacobian','gradxB','gradyB','gradzB','gradz_coeff',...
         'gxx','gxy','gyy','gyz','gzz','hatB','hatR','hatZ'};
geo_arrays = zeros(2,data.Nz,numel(names));

for i_ = 1:numel(names)
    namae = names{i_};
    geo_arrays(:,:,i_) = h5read(data.outfilenames{end},['/data/metric/',namae])';
end
fig = figure;
subplot(311)
    for i = 1:5
    plot(data.z, geo_arrays(1,:,i),'DisplayName',names{i}); hold on;
    end
    xlim([min(data.z),max(data.z)]); legend('show'); title('MoNoLiT geometry');

subplot(312)
    for i = 6:10
    plot(data.z, geo_arrays(1,:,i),'DisplayName',names{i}); hold on;
    end
    xlim([min(data.z),max(data.z)]); legend('show');

subplot(313)
    for i = 11:13
    plot(data.z, geo_arrays(1,:,i),'DisplayName',names{i}); hold on;
    end
    xlim([min(data.z),max(data.z)]); legend('show');
end

