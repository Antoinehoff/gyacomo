function [ FIGURE ] = spectrum_1D( data, options )

options.PLAN      = 'kxky';
options.COMP      = 'avg';
options.INTERP    = 0;
options.POLARPLOT = 0;
options.AXISEQUAL = 1;
toplot = process_field(data,options);
t = data.Ts3D; frames = toplot.FRAMES;

colors = jet(numel(frames));

switch options.NAME
    case '\Gamma_x'
    FIGURE.fig = figure; FIGURE.FIGNAME = ['transp_spectrum_',data.PARAMS]; 
    yname = '$\sum_{k_y}|\Gamma_k|$'; fieldname = 'transport';
    case '\phi'
    FIGURE.fig = figure; FIGURE.FIGNAME = ['phi_spectrum_',data.PARAMS]; 
    yname = '$\sum_{k_y}|\phi|$'; fieldname = 'ES pot.';
    case 'n_i'
    FIGURE.fig = figure; FIGURE.FIGNAME = ['ni_spectrum_',data.PARAMS]; 
    yname = '$\sum_{k_y}|n_i|$'; fieldname = 'ion dens.';
    case 'n_e'
    FIGURE.fig = figure; FIGURE.FIGNAME = ['ne_spectrum_',data.PARAMS]; 
    yname = '$\sum_{k_y}|n_e|$'; fieldname = 'elec. dens.';
end

PLOT2D = 0;
switch options.COMPXY
    case 'avg'
        compx = @(x) mean(x,2);
        compy = @(x) mean(x,1);
    case 'sum'
        compx = @(x) sum(x,2);
        compy = @(x) sum(x,1);
    case 'max'
        compx = @(x) max(x,2);
        compy = @(x) max(x,1);
    otherwise
        compx =  @(x) x(:,:);
        compy =  @(x) x(:,:);
        PLOT2D= 1;
end

if ~PLOT2D
    set(gcf, 'Position',  [20 50 1200 500])
    subplot(1,2,1)

    k     = data.ky;
    xname = '$k_y$';

    nmax  = ceil(data.Nkx*2/3);
    shiftx = @(x) x;%(1:nmax);
    shifty = @(x) x;%(1:nmax);
    switch options.COMPT
    case 'avg'
        it0 = toplot.FRAMES(1); it1 = toplot.FRAMES(end);
        Gk    = compx(abs(mean(toplot.FIELD(:,:,:),3)));
        Gk    = squeeze(Gk);
        if options.NORM
            Gk    = Gk./max(max(abs(Gk)));
        end
        X = shiftx(k);
        if options.OK
            Y = shifty(Gk)./X;
        else
            Y = shifty(Gk);
        end
        plot(X,Y,'DisplayName','t-averaged')  
    otherwise
    for it = 1:numel(frames)
        Gk    = compx(abs(toplot.FIELD(:,:,it)));
        Gk    = squeeze(Gk);
        if options.NORM
            Gk    = Gk./max(max(abs(Gk)));
        end
        X = shiftx(k);
        if options.OK
            Y = shifty(Gk)./X;
        else
            Y = shifty(Gk);
        end
        plot(X,Y,'DisplayName',['$t=',num2str(t(frames(it))),'$'],...
            'Color',colors(it,:)); hold on;
    end
    end
    grid on
    title(['HeLaZ $k_x$ ',fieldname,' spectrum']); legend('show','Location','eastoutside')
    xlabel(xname); ylabel(yname)

    subplot(1,2,2)

    k     = data.kx;
    xname = '$k_x$';
    nmax  = floor(data.Nky/2*2/3);
    switch options.COMPT
    case 'avg'
%         it0 = toplot.FRAMES(1); it1 = toplot.FRAMES(end);
        Gk    = compy(abs(mean(toplot.FIELD(:,:,:),3)))';
        Gk    = squeeze(Gk);
        if options.NORM
            Gk    = Gk./max(max(abs(Gk)));
        end
        X = k(k>0); X = X(1:end-1);
        Yp= Gk(k>0);
        Ym= Gk(k<0);
        Y = Yp(1:end-1) + Ym(end:-1:1); 
        Y = Y(end:-1:1);
        plot(X,Y,'DisplayName','t-averaged')  
    otherwise
    for it = 1:numel(toplot.FRAMES)
        Gk    = compy(abs(toplot.FIELD(:,:,it)));
        Gk    = squeeze(Gk);
        if options.NORM
            Gk    = Gk./max(max(abs(Gk)));
        end
        X = k(k>0); X = X(1:end-1);
        Yp= Gk(k>0);
        Ym= Gk(k<0);
        Y = Yp(1:end-1) + Ym(end:-1:1); 
        Y = Y(end:-1:1);
        plot(X,Y,'DisplayName',['$t=',num2str(t(frames(it))),'$'],...
            'Color',colors(it,:)); hold on;
    end
    end
    grid on
    %     title('HeLaZ $k_y$ transport spectrum'); 
    legend('show','Location','eastoutside');
    xlabel(xname); ylabel(yname)
else
%     it0 = toplot.FRAMES(1); it1 = toplot.FRAMES(end);
    Gk    = mean(abs(toplot.FIELD(:,:,:)),3);
    Gk    = squeeze(Gk);
    if options.NORM
        Gk    = Gk./max(max(abs(Gk)));
    end
    pclr = pcolor(toplot.X,toplot.Y,Gk);
    set(pclr, 'edgecolor','none');

end

