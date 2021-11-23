function [ FIGURE ] = photomaton( DATA,OPTIONS )
%UNTITLED5 Summary of this function goes here
%% Processing
toplot = process_field(DATA,OPTIONS);
FNAME  = toplot.FILENAME;
FRAMES = toplot.FRAMES;

%
TNAME = [];
FIGURE.fig = figure; set(gcf, 'Position',  toplot.DIMENSIONS.*[1 1 numel(FRAMES) 1])
    for i_ = 1:numel(FRAMES)
    subplot(1,numel(FRAMES),i_); TNAME = [TNAME,'_',sprintf('%.0f',DATA.Ts3D(FRAMES(i_)))];
        if ~strcmp(OPTIONS.PLAN,'kxky')
            scale = max(max(abs(toplot.FIELD(:,:,FRAMES(i_))))); % Scaling to normalize
        else
            scale = 1;
        end
        pclr = pcolor(toplot.X,toplot.Y,toplot.FIELD(:,:,FRAMES(i_))./scale); set(pclr, 'edgecolor','none');pbaspect(toplot.ASPECT)
%         if ~strcmp(OPTIONS.PLAN,'kxky')
%             caxis([-1,1]);
%             colormap(bluewhitered);
%         end
        xlabel(toplot.XNAME); ylabel(toplot.YNAME);
%         if i_ > 1; set(gca,'ytick',[]); end; 
        title([sprintf('$t c_s/R=%.0f$',DATA.Ts3D(FRAMES(i_))),', max = ',sprintf('%.1e',scale)]);
        if OPTIONS.INTERP
            shading interp; 
        end
    end
    legend(['$',toplot.FIELDNAME,'$']);
    FIGURE.FIGNAME = [FNAME,'_snaps',TNAME];
end

