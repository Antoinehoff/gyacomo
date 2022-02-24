function [ FIGURE ] = real_plot_1D( data, options )

options.PLAN      = 'xy';
options.COMP      = options.COMPZ;
options.POLARPLOT = 0;
options.INTERP    = 0;

toplot = process_field(data,options);
t = data.Ts3D; frames = toplot.FRAMES;

colors = jet(numel(frames));

switch options.NAME
    case '\Gamma_x'
    FIGURE.fig = figure; FIGURE.FIGNAME = ['transport_1D_avg_',data.PARAMS]; 
    yname = '\Gamma';
    case 'v_y'
    FIGURE.fig = figure; FIGURE.FIGNAME = ['ZF_1D_avg_',data.PARAMS]; 
    yname = 'v_y';
    case 's_y'
    FIGURE.fig = figure; FIGURE.FIGNAME = ['ZS_1D_avg_',data.PARAMS]; 
    yname = 's_y';
    case '\phi'
    FIGURE.fig = figure; FIGURE.FIGNAME = ['phi_1D_avg_',data.PARAMS]; 
    yname = '\phi';
    case 'n_i'
    FIGURE.fig = figure; FIGURE.FIGNAME = ['ni_1D_avg_',data.PARAMS]; 
    yname = 'n_i';
    case 'n_e'
    FIGURE.fig = figure; FIGURE.FIGNAME = ['ne_1D_avg_',data.PARAMS]; 
    yname = 'n_e';
end

switch options.COMPXY
    case 'avg'
        compx  = @(x) mean(x,1);
        compy  = @(x) mean(x,2);
        ynamex = ['$\langle ',yname,'\rangle_x$'];
        ynamey = ['$\langle ',yname,'\rangle_y$'];
    otherwise
        compx  = @(x) x(1,:);
        compy  = @(x) x(:,1);
        ynamex = ['$',yname,'(x=0)$'];
        ynamey = ['$',yname,'(y=0)$'];
end
    

set(gcf, 'Position',  [20 50 5000 2000])
subplot(1,2,1)

    X     = data.x;
    xname = '$x$';
    switch options.COMPT
        case 'avg'
            Y    = compy(mean(toplot.FIELD(:,:,:),3));
            Y    = squeeze(Y);
            if options.NORM
                Y    = Y./max(abs(Y));
            end    
            plot(X,Y,'DisplayName',['$t\in[',num2str(t(frames(1))),',',num2str(t(frames(end))),']$'],...
            'Color',colors(1,:)); hold on;
            legend('show')
        otherwise
            for it = 1:numel(toplot.FRAMES)
                Y    = compy(toplot.FIELD(:,:,toplot.FRAMES(it)));
                Y    = squeeze(Y);
                if options.NORM
                    Y    = Y./max(abs(Y));
                end
                plot(X,Y,'DisplayName',['$t=',num2str(t(frames(it))),'$'],...
                    'Color',colors(it,:)); hold on;
            end
            legend('show','Location','eastoutside')
    end
    grid on
    xlabel(xname); ylabel(ynamey)
    xlim([min(X),max(X)]);
    
subplot(1,2,2)

    X     = data.y;
    xname = '$y$';
    switch options.COMPT
        case 'avg'
            Y    = compx(mean(toplot.FIELD(:,:,:),3));
            Y    = squeeze(Y);
            if options.NORM
                Y    = Y./max(abs(Y));
            end    
            plot(X,Y,'DisplayName',['$t\in[',num2str(t(frames(1))),',',num2str(t(frames(end))),']$'],...
            'Color',colors(1,:)); hold on;
            legend('show')
        otherwise
            for it = 1:numel(toplot.FRAMES)
                Y    = compx(toplot.FIELD(:,:,toplot.FRAMES(it)));
                Y    = squeeze(Y);
                if options.NORM
                    Y    = Y./max(abs(Y));
                end
                plot(X,Y,'DisplayName',['$t=',num2str(t(frames(it))),'$'],...
                    'Color',colors(it,:)); hold on;
            end
            legend('show','Location','eastoutside')
    end

    grid on
%     title('HeLaZ $k_y$ transport spectrum'); 
    legend('show','Location','eastoutside');
    xlabel(xname); ylabel(ynamex)
        xlim([min(X),max(X)]);

end

