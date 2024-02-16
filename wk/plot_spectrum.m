function [FIGURE] = plot_spectrum(data,options)


options.INTERP    = 0;
options.POLARPLOT = 0;
options.AXISEQUAL = 1;
options.LOGSCALE  = 0;
options.CLIMAUTO  = 1;
options.PLAN      = 'kxky'; options.COMP = data.grids.Nz/2+1;
options.RESOLUTION = 256;
options.BWR       = 0; % bluewhitered plot or gray
options.COMP      = 'avg';
[toplot] = process_field(data,options);

% Time averaging
fkykx = squeeze(mean(abs(toplot.FIELD(:,:,:)),3));

FNAME  = toplot.FILENAME;


kx  = toplot.X(1,:);
Nkx = numel(kx);
kx  = toplot.X(1,Nkx/2:end);
fkx = squeeze(fkykx(1,Nkx/2:end));
ky  = toplot.Y(:,data.grids.Nkx/2+1);
fky = squeeze(fkykx(:,data.grids.Nkx/2+1));

FIGURE.fig = figure;
subplot(121)
    scale = max(fkx)*options.NORMALIZE + 1*(1-options.NORMALIZE);
    plot(kx,fkx./scale)
    xlabel(toplot.XNAME)
    ylabel(['$|',toplot.FIELDNAME,'|$'])
    

subplot(122)
    scale = max(fky)*options.NORMALIZE + 1*(1-options.NORMALIZE);
    plot(ky,fky./scale)
    xlabel(toplot.YNAME)
    ylabel(['$|',toplot.FIELDNAME,'|$'])


FIGURE.FIGNAME = [FNAME,'_spectrum'];

end