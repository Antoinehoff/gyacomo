title1 = GIFNAME;
FIGDIR = BASIC.RESDIR;
if ~exist(FIGDIR, 'dir')
   mkdir(FIGDIR)
end

GIFNAME = [FIGDIR, GIFNAME,'.gif'];

% Set colormap boundaries
hmax = max(max(max(FIELD)));
hmin = min(min(min(FIELD)));
flag = 0;
if hmax == hmin 
    disp('Warning : h = hmin = hmax = const')
else
% Setup figure frame
fig  = figure;
    plot(X,FIELD(:,1)); % to set up
    axis tight manual % this ensures that getframe() returns a consistent size
    in      = 1;
    nbytes = fprintf(2,'frame %d/%d',in,numel(FIELD(1,1,:)));
    for n = FRAMES % loop over selected frames
        scale = max(FIELD(:,n))*SCALING + (1-SCALING);
        plot(X,FIELD(:,n)/scale,linestyle);
        if (YMIN ~= YMAX && XMIN ~= XMAX)
        ylim([YMIN,YMAX]); xlim([XMIN,XMAX]);
        end
        title(['$t \approx$', sprintf('%.3d',ceil(T(n))), ', scaling = ',sprintf('%.1e',scale)]);
        xlabel(XNAME); ylabel(FIELDNAME);
        drawnow 
        % Capture the plot as an image 
        frame = getframe(fig); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,64); 
        % Write to the GIF File 
        if in == 1 
          imwrite(imind,cm,GIFNAME,'gif', 'Loopcount',inf); 
        else 
          imwrite(imind,cm,GIFNAME,'gif','WriteMode','append', 'DelayTime',DELAY); 
        end 
        % terminal info
        while nbytes > 0
          fprintf('\b')
          nbytes = nbytes - 1;
        end
        nbytes = fprintf(2,'frame %d/%d',n,numel(FIELD(1,1,:)));
        in = in + 1;
    end
    disp(' ')
    disp(['Gif saved @ : ',GIFNAME])
end

