function [ FIGURE ] = photomaton( DATA,OPTIONS )
%UNTITLED5 Summary of this function goes here
%% Processing
switch OPTIONS.PLAN
    case 'RZ'
        toplot = poloidal_plot(DATA,OPTIONS);
    otherwise
        toplot = process_field(DATA,OPTIONS);
end
FNAME  = toplot.FILENAME;
FRAMES = toplot.FRAMES;
Nframes= numel(FRAMES);
Nrows  = ceil(Nframes/4);
Ncols  = ceil(Nframes/Nrows);
%
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
            pclr = pcolor(toplot.X,toplot.Y,toplot.FIELD(:,:,i_)./scale); set(pclr, 'edgecolor','none');
            if OPTIONS.AXISEQUAL
                pbaspect(toplot.ASPECT)
            end
            if ~strcmp(OPTIONS.PLAN,'kxky')
                clim([-1,1]*frame_max/scale);
                colormap(bluewhitered);
                if OPTIONS.INTERP
                    shading interp; 
                end
            end
        else
            contour(toplot.X,toplot.Y,toplot.FIELD(:,:,i_)./scale,128);
%             pclr = pcolor(toplot.X,toplot.Y,toplot.FIELD(:,:,FRAMES(i_))./scale); set(pclr, 'edgecolor','none'); shading interp
            tshot = DATA.Ts5D(FRAMES(i_));
        end
        xlabel(toplot.XNAME); ylabel(toplot.YNAME);
%         if i_ > 1; set(gca,'ytick',[]); end; 
        title([sprintf('$t c_s/R=%.0f$',tshot),', max = ',sprintf('%.1e',frame_max)]);
        if OPTIONS.CLIMAUTO
            clim('auto')
        end
    end
    legend(['$',toplot.FIELDNAME,'$']);
    FIGURE.FIGNAME = [FNAME,'_snaps',TNAME];
end

