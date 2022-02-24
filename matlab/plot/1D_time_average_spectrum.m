function [ FIGURE ] = spectrum_1D( data, options )

options.PLAN      = 'kxky';
options.COMP      = 1;

toplot = process_field(data,options);
t = data.Ts3D; frames = toplot.FRAMES;

colors = jet(numel(frames));

switch options.NAME
    case '\Gamma_x';
    FIGURE.fig = figure; FIGURE.FIGNAME = ['transp_spectrum_',data.PARAMS]; 
    case '\phi'
    FIGURE.fig = figure; FIGURE.FIGNAME = ['phi_spectrum_',data.PARAMS]; 
    case 'n_i'
    FIGURE.fig = figure; FIGURE.FIGNAME = ['ni_spectrum_',data.PARAMS]; 
    case 'n_e'
    FIGURE.fig = figure; FIGURE.FIGNAME = ['ne_spectrum_',data.PARAMS]; 
end
set(gcf, 'Position',  [20 50 5000 2000])
subplot(1,2,1)

    sdim  = 2;
    k     = data.kx;
    xname = '$k_x$';
    yname = '$\sum_{k_y}|\Gamma_k|$';
    nmax  = ceil(data.Nkx*2/3);
    shiftx = @(x) x(1:nmax);
    shifty = @(x) x(1:nmax);

    for it = 1:numel(toplot.FRAMES)
        Gk    = sum(abs(toplot.FIELD(:,:,toplot.FRAMES(it))),sdim);
        Gk    = squeeze(Gk);
        if options.NORM
            Gk    = Gk./max(abs(Gk));
        end
        semilogy(shiftx(k), shifty(Gk),'DisplayName',['$t=',num2str(t(frames(it))),'$'],...
            'Color',colors(it,:)); hold on;
    end
    grid on
    title('HeLaZ $k_x$ transport spectrum'); legend('show','Location','eastoutside')
    xlabel(xname); ylabel(yname)
    
subplot(1,2,2)

    sdim  = 1;
    k     = data.ky;
    xname = '$k_y$';
    yname = '$\sum_{k_x}|\Gamma_k|$';
    nmax  = floor(data.Nky/2*2/3);
    AA     = @(x) x(1:nmax);
    shiftx = @(x) AA(x);
    shifty = @(x) AA(ifftshift(x));

    for it = 1:numel(toplot.FRAMES)
        Gk    = sum(abs(toplot.FIELD(:,:,toplot.FRAMES(it))),sdim);
        Gk    = squeeze(Gk);
        if options.NORM
            Gk    = Gk./max(abs(Gk));
        end
        semilogy(shiftx(k), shifty(Gk),'DisplayName',['$t=',num2str(t(frames(it))),'$'],...
            'Color',colors(it,:)); hold on;
    end
    grid on
    title('HeLaZ $k_y$ transport spectrum'); legend('show','Location','eastoutside');
    xlabel(xname); ylabel(yname)
end

