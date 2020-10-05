title1 = GIFNAME;
FIGDIR = ['../results/',BASIC.SIMID,'/'];
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
    pcolor(X,Y,FIELD(:,:,1)); % to set up
    colormap jet
    colorbar
    axis tight manual % this ensures that getframe() returns a consistent size
    shading interp;
    xlabel(XNAME); ylabel(YNAME);
    hold on
  
    in      = 1;
    nbytes = fprintf(2,'frame %d/%d',in,numel(FIELD(1,1,:)));
    for n = FRAMES % loop over selected frames
        pclr = pcolor(X,Y,FIELD(:,:,n)); shading interp; % frame plot
        set(pclr, 'edgecolor','none');
        title([FIELDNAME,', $t \approx$', sprintf('%.3d',ceil(T(n)))]);
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

