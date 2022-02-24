function [ FIGURE ] = spectrum_1D( data, options )

options.PLAN      = 'kxky';
options.COMP      = 1;
options.POLARPLOT = 0;
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

switch options.COMPXY
    case 'avg'
        compx = @(x) mean(x,1);
        compy = @(x) mean(x,2);
    case 'sum'
        compx = @(x) sum(x,1);
        compy = @(x) sum(x,2);
    case 'max'
        compx = @(x) max(x,1);
        compy = @(x) max(x,2);
    otherwise
        compx =  @(x) x(data.kx==0,:);
        compy =  @(x) x(:,data.ky==0);
end

set(gcf, 'Position',  [20 50 5000 2000])
subplot(1,2,1)

    k     = data.kx;
    xname = '$k_x$';

    nmax  = ceil(data.Nkx*2/3);
    shiftx = @(x) x;%(1:nmax);
    shifty = @(x) x;%(1:nmax);
    switch options.COMPT
        case 'avg'
            it0 = toplot.FRAMES(1); it1 = toplot.FRAMES(end);
            Gk    = compy(abs(mean(toplot.FIELD(:,:,it0:it1),3)));
            Gk    = squeeze(Gk);
            if options.NORM
                Gk    = Gk./max(abs(Gk));
            end
            X = shiftx(k);
            if options.OK
                Y = shifty(Gk)./X;
            else
                Y = shifty(Gk);
            end
            plot(X,Y,'DisplayName','t-averaged')  
        otherwise
        for it = 1:numel(toplot.FRAMES)
            Gk    = compy(abs(toplot.FIELD(:,:,toplot.FRAMES(it))));
            Gk    = squeeze(Gk);
            if options.NORM
                Gk    = Gk./max(abs(Gk));
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

    k     = data.ky;
    xname = '$k_y$';
    nmax  = floor(data.Nky/2*2/3);
    AA     = @(x) x(1:nmax);
    shiftx = @(x) x;%AA(x);
    shifty = @(x) x;%AA(ifftshift(x));
    switch options.COMPT
        case 'avg'
            it0 = toplot.FRAMES(1); it1 = toplot.FRAMES(end);
            Gk    = compx(abs(mean(toplot.FIELD(:,:,it0:it1),3)))';
            Gk    = squeeze(Gk);
            if options.NORM
                Gk    = Gk./max(abs(Gk));
            end
            X = k(k>0); X = X(1:end-1);
            Yp= Gk(k>0);
            Ym= Gk(k<0);
            Y = Yp(1:end-1) + Ym(end:-1:1); 
            Y = Y(end:-1:1);
            plot(X,Y,'DisplayName','t-averaged')  
        otherwise
        for it = 1:numel(toplot.FRAMES)
            Gk    = compx(abs(toplot.FIELD(:,:,toplot.FRAMES(it))));
            Gk    = squeeze(Gk);
            if options.NORM
                Gk    = Gk./max(abs(Gk));
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
end

