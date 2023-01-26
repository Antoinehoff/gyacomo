function [full_fig] = dashboard(DATA)
% Make an overview of the simulations results in one figure
full_fig = figure;
axes = 1:9;

axes(1) = subplot(8,1,1,'parent',full_fig);
axes(2) = subplot(8,1,2,'parent',full_fig);
axes(3) = subplot(4,2,3,'parent',full_fig);
axes(4) = subplot(4,4,7,'parent',full_fig);
axes(5) = subplot(4,4,8,'parent',full_fig);
axes(6) = subplot(4,2,5,'parent',full_fig);
axes(7) = subplot(4,2,6,'parent',full_fig);
axes(8) = subplot(4,2,7,'parent',full_fig);
axes(9) = subplot(4,2,8,'parent',full_fig);


%% Space time diagramm (fig 11 Ivanov 2020)
options.TAVG_0   = 0.1*DATA.Ts3D(end);
options.TAVG_1   = DATA.Ts3D(end); % Averaging times duration
options.NMVA     = 1;              % Moving average for time traces
options.ST_FIELD = '\phi';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x, Q_x)
options.INTERP   = 1;
options.NCUT     = 4;              % Number of cuts for averaging and error estimation
options.RESOLUTION = 256;
PLOT = plot_radial_transport_and_spacetime(DATA,options,'GENE');
% Put on full fig
axcp = copyobj(PLOT.ax1,full_fig); 
set(axcp,'Position',get(axes(1),'position'));delete(axes(1));
axcp = copyobj(PLOT.ax3,full_fig); 
set(axcp,'Position',get(axes(2),'position'));delete(axes(2));
colormap(axcp,bluewhitered);
close(PLOT.fig);

%% Show f_i(vpar,mu)
options.T         = [ 1]*DATA.Ts3D(end);
options.SPECIES   = 'i';
% options.PLT_FCT = 'contour';
options.PLT_FCT = 'contourf';
% options.PLT_FCT = 'surf';
% options.PLT_FCT = 'surfvv';
options.non_adiab = 0;
options.RMS       = 1; % Root mean square i.e. sqrt(sum_k|f_k|^2) as in Gene
options.folder  = DATA.folder;
options.iz      = 'avg';
options.FIELD   = '<f_>';
options.SPAR    = linspace(-3,3,32);
options.XPERP   = linspace( 0,sqrt(6),16).^2;

options.ONED    = 0;
switch DATA.CODENAME
    case 'GENE'
        PLOT = plot_fa_gene(options);
    case 'GYACOMO'
        PLOT = plot_fa(DATA,options);
end
% Put on full fig
axcp = copyobj(PLOT.ax1,full_fig); 
set(axcp,'Position',get(axes(3),'position'));delete(axes(3));
close(PLOT.fig);

options.ONED    = 1;
switch DATA.CODENAME
    case 'GENE'
        PLOT = plot_fa_gene(options);
    case 'GYACOMO'
        PLOT = plot_fa(DATA,options);
end
% Put on full fig
axcp = copyobj(PLOT.ax1,full_fig); 
set(axcp,'Position',get(axes(4),'position'));delete(axes(4));
axcp = copyobj(PLOT.ax2,full_fig); 
set(axcp,'Position',get(axes(5),'position'));delete(axes(5));
close(PLOT.fig);

%% Time averaged spectrum
options.TIME   = [100 500];
options.NORM   =1;
% options.NAME   = '\phi';
% options.NAME      = 'n_i';
options.NAME   ='\Gamma_x';
options.PLAN   = 'kxky';
options.COMPZ  = 'avg';
options.OK     = 0;
options.COMPXY = 'avg'; % avg/sum/max/zero/ 2D plot otherwise
options.COMPT  = 'avg';
options.PLOT   = 'semilogy';
PLOT = spectrum_1D(DATA,options);

% Put on full fig
axcp = copyobj(PLOT.ax1,full_fig); 
set(axcp,'Position',get(axes(6),'position'));delete(axes(6));
axcp = copyobj(PLOT.ax2,full_fig); 
set(axcp,'Position',get(axes(7),'position'));delete(axes(7));
close(PLOT.fig);

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
PLOT = mode_growth_meter(DATA,options);

% Put on full fig
axcp = copyobj(PLOT.axes(1),full_fig); 
set(axcp,'Position',get(axes(8),'position'));delete(axes(8));
axcp = copyobj(PLOT.axes(4),full_fig); 
set(axcp,'Position',get(axes(9),'position'));delete(axes(9));
set(gca,'xtick',[]);
close(PLOT.fig);
end
