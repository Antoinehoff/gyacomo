title1 = GIFNAME;
FIGDIR = BASIC.RESDIR;
GIFNAME = [FIGDIR, GIFNAME];
XNAME     = latexize(XNAME);
YNAME     = latexize(YNAME);
FIELDNAME = latexize(FIELDNAME);

vidfile = VideoWriter(GIFNAME,'Uncompressed AVI');
vidfile.FrameRate = FPS;
open(vidfile);
% Set colormap boundaries
hmax = max(max(max(FIELD(:,:,FRAMES))));
hmin = min(min(min(FIELD(:,:,FRAMES))));
scale = -1;
flag = 0;
if hmax == hmin 
    disp('Warning : h = hmin = hmax = const')
else
% Setup figure frame
figure('Color','white','Position', [100, 100, 500, 500]);
    pcolor(X,Y,FIELD(:,:,1)); % to set up
%     colormap gray
    colormap(bluewhitered)
    axis tight manual % this ensures that getframe() returns a consistent size
    if INTERP
        shading interp;
    end
    in      = 1;
    nbytes = fprintf(2,'frame %d/%d',in,numel(FIELD(1,1,:)));
    for n = FRAMES % loop over selected frames
        scale = max(max(abs(FIELD(:,:,n))));
        pclr = pcolor(X,Y,FIELD(:,:,n)/scale); % frame plot
        if INTERP
            shading interp; 
        end
        if NORMALIZED
            pclr = pcolor(X,Y,FIELD(:,:,n)/scale); % frame plot\
            caxis([-1,1]);
        else
            pclr = pcolor(X,Y,FIELD(:,:,n)); % frame plot\
            if CONST_CMAP
                caxis([-1,1]*max(abs([hmin hmax]))); % const color map
            else
                caxis([-1,1]*scale); % adaptive color map                
            end
            title([FIELDNAME,', $t \approx$', sprintf('%.3d',ceil(T(n)))...
                ,', scale = ',sprintf('%.1e',scale)]);
        end
        set(pclr, 'edgecolor','none'); axis square;
        xlabel(XNAME); ylabel(YNAME); %colorbar;
        title([FIELDNAME,', $t \approx$', sprintf('%.3d',ceil(T(n)))...
            ,', scaling = ',sprintf('%.1e',scale)]);
        drawnow 
        % Capture the plot as an image 
        frame = getframe(gcf); 
        writeVideo(vidfile,frame);
        % terminal info
        while nbytes > 0
          fprintf('\b')
          nbytes = nbytes - 1;
        end
        nbytes = fprintf(2,'frame %d/%d',n,numel(FIELD(1,1,:)));
        in = in + 1;
    end
    disp(' ')
    disp(['Video saved @ : ',GIFNAME])
    close(vidfile);
end

