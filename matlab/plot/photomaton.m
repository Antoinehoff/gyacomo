function [ FIGURE ] = photomaton( DATA,OPTIONS )
%UNTITLED5 Summary of this function goes here
%% Processing
switch OPTIONS.PLAN
    case 'RZ'
        toplot = poloidal_plot(DATA,OPTIONS);
    case '3D'
        OPTIONS.PLAN      = 'xy';  
        [~,OPTIONS.COMP] = min(abs(OPTIONS.XYZ(3) - DATA.grids.z));
        toplot    = process_field(DATA,OPTIONS);
        OPTIONS.PLAN      = 'xz';  
        [~,OPTIONS.COMP] = min(abs(OPTIONS.XYZ(2) - DATA.grids.y));
        toplot_xz = process_field(DATA,OPTIONS);
        OPTIONS.PLAN      = 'yz';  
        [~,OPTIONS.COMP] = min(abs(OPTIONS.XYZ(1) - DATA.grids.x));
        toplot_yz = process_field(DATA,OPTIONS);
        OPTIONS.PLAN      = '3D';  
    otherwise
        toplot = process_field(DATA,OPTIONS);
end
FNAME  = toplot.FILENAME;
FRAMES = toplot.FRAMES;
if OPTIONS.TAVG
    toplot.FIELD = mean(toplot.FIELD,3);
    FRAMES = FRAMES(1);
end
Nframes= numel(FRAMES);
Nrows  = ceil(Nframes/4);
Ncols  = ceil(Nframes/Nrows);
%
if OPTIONS.LOGSCALE
    toplot.FIELD = max(abs(toplot.FIELD),1e-40);
    toplot.FIELD = log(toplot.FIELD);
end
TNAME = [];
FIGURE.fig = figure; %set(gcf, 'Position',  toplot.DIMENSIONS.*[1 1 Ncols Nrows])
    for i_ = 1:numel(FRAMES)
    frame_max = max(max(max(abs(toplot.FIELD(:,:,i_)))));
    subplot(Nrows,Ncols,i_); TNAME = [TNAME,'_',sprintf('%.0f',DATA.Ts3D(FRAMES(i_)))];
    if OPTIONS.NORMALIZE
        scale = frame_max; % Scaling to normalize
    else
        scale = 1;
    end
        if ~strcmp(OPTIONS.PLAN,'sx')
            tshot = DATA.Ts3D(FRAMES(i_));
            if ~strcmp(OPTIONS.PLAN,'3D')
                pclr = pcolor(toplot.X,toplot.Y,toplot.FIELD(:,:,i_)./scale); set(pclr, 'edgecolor','none');
            else
                % xy plane
               s = surface(toplot.X,toplot.Y,OPTIONS.XYZ(3)+0*toplot.Y,toplot.FIELD(:,:,i_)./scale); hold on;
                s.EdgeColor = 'none';
                % s.FaceAlpha = 0.5;
                % xz plane
               s = surface(toplot_xz.X,OPTIONS.XYZ(2)+0*toplot_xz.X,toplot_xz.Y,toplot_xz.FIELD(:,:,i_)./scale);
                s.EdgeColor = 'none';
                % s.FaceAlpha = 0.7;
                % yz plane
               s = surface(OPTIONS.XYZ(1)+0*toplot_yz.X,toplot_yz.X,toplot_yz.Y,toplot_yz.FIELD(:,:,i_)./scale);
                s.EdgeColor = 'none';
                % s.FaceAlpha = 0.7;
               xlabel('x');
               ylabel('y');
               zlabel('z');
               h = gca;
               plot3(h.XLim,OPTIONS.XYZ(2)*[1 1],OPTIONS.XYZ(3)*[1 1],'--k','LineWidth',1.)
               plot3(OPTIONS.XYZ(1)*[1 1],h.YLim,OPTIONS.XYZ(3)*[1 1],'--k','LineWidth',1.)
               plot3(OPTIONS.XYZ(1)*[1 1],OPTIONS.XYZ(2)*[1 1],h.ZLim,'--k','LineWidth',1.)
               % zline(OPTIONS.XYZ(3),'-k','LineWidth',1.5)
               view([1 -1 0.25])
           end
        else
            contour(toplot.X,toplot.Y,toplot.FIELD(:,:,i_)./scale,128);
%             pclr = pcolor(toplot.X,toplot.Y,toplot.FIELD(:,:,FRAMES(i_))./scale); set(pclr, 'edgecolor','none'); shading interp
            tshot = DATA.Ts5D(FRAMES(i_));
        end
        xlabel(toplot.XNAME); ylabel(toplot.YNAME);
        title([sprintf('$t c_s/R=%5.2f$',tshot),', max = ',sprintf('%.1e',frame_max)]);
    end
    if OPTIONS.AXISEQUAL
        pbaspect(toplot.ASPECT)
    end
    if ~strcmp(OPTIONS.PLAN,'kxky')
        % clim([-1,1]*frame_max/scale);
        colormap(bluewhitered);
        if OPTIONS.INTERP
            shading interp; 
        end
    end
    legend(['$',toplot.FIELDNAME,'$']);
    legend('Location','northeast');
    FIGURE.FIGNAME = [FNAME,'_snaps',TNAME];
end

