title1 = GIFNAME;
FIGDIR = BASIC.RESDIR;
GIFNAME = [FIGDIR, GIFNAME,'.mp4'];
vidfile = VideoWriter(GIFNAME,'Uncompressed AVI');
open(vidfile);
% Set colormap boundaries
hmax = max(max(max(FIELD)));
hmin = min(min(min(FIELD)));

flag = 0;
if hmax == hmin 
    disp('Warning : h = hmin = hmax = const')
else
% Setup figure frame
figure('Color','white','Position', [100, 100, 400, 400]);
    pcolor(X,Y,FIELD(:,:,1)); % to set up
    colormap gray
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
        set(pclr, 'edgecolor','none'); axis square;
%         caxis([-max(max(abs(FIELD(:,:,n)))),max(max(abs(FIELD(:,:,n))))]);
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

