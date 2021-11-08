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
        pclr = pcolor(toplot.X,toplot.Y,toplot.FIELD(:,:,FRAMES(i_))); set(pclr, 'edgecolor','none');pbaspect(toplot.ASPECT)
        colormap(bluewhitered);
        xlabel(toplot.XNAME); ylabel(toplot.YNAME);set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',DATA.Ts3D(FRAMES(i_))));
        if OPTIONS.INTERP
            shading interp; 
        end
    end
    legend(['$',toplot.FIELDNAME,'$']);
    FIGURE.FIGNAME = [FNAME,'_snaps',TNAME];
end

