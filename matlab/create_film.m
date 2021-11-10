function create_film(DATA,OPTIONS,format)
%% Plot options
FPS = 30; DELAY = 1/FPS;
BWR = 1; NORMALIZED = 1;
T = DATA.Ts3D;
%% Processing
toplot = process_field(DATA,OPTIONS);

%%
FILENAME  = [DATA.localdir,toplot.FILENAME,format];
XNAME     = ['$',toplot.XNAME,'$'];
YNAME     = ['$',toplot.YNAME,'$'];
FIELDNAME = ['$',toplot.FIELDNAME,'$'];
FIELD     = toplot.FIELD; X = toplot.X; Y = toplot.Y;
FRAMES    = toplot.FRAMES;
switch format
    case '.avi'
        vidfile = VideoWriter(FILENAME,'Uncompressed AVI');
        vidfile.FrameRate = FPS;
        open(vidfile);  
end
% Set colormap boundaries
hmax = max(max(max(FIELD(:,:,FRAMES))));
hmin = min(min(min(FIELD(:,:,FRAMES))));
scale = -1;
flag = 0;
if hmax == hmin 
    disp('Warning : h = hmin = hmax = const')
else
% Setup figure frame
fig  = figure('Color','white','Position', toplot.DIMENSIONS.*[1 1 1.2 1]);
    pcolor(X,Y,FIELD(:,:,1)); % to set up
    if BWR
    colormap(bluewhitered)
    end
    axis tight manual % this ensures that getframe() returns a consistent size
    if toplot.INTERP
        shading interp;
    end
    in      = 1;
    nbytes = fprintf(2,'frame %d/%d',in,numel(FIELD(1,1,:)));
    for n = FRAMES % loop over selected frames
        scale = max(max(abs(FIELD(:,:,n)))); % Scaling to normalize
        if NORMALIZED
            pclr = pcolor(X,Y,FIELD(:,:,n)/scale); % frame plot\
            caxis([-1,1]);
        else
            pclr = pcolor(X,Y,FIELD(:,:,n)); % frame plot\
            if CONST_CMAP
                caxis([-1,1]*max(abs([hmin hmax]))); % adaptive color map
            else
                caxis([-1,1]*scale); % adaptive color map                
            end
            title([FIELDNAME,', $t \approx$', sprintf('%.3d',ceil(T(n)))...
                ,', scale = ',sprintf('%.1e',scale)]);
        end
        title([FIELDNAME,', $t \approx$', sprintf('%.3d',ceil(T(n)))...
              ,', scaling = ',sprintf('%.1e',scale)]);
        if toplot.INTERP
            shading interp; 
        end
        set(pclr, 'edgecolor','none'); pbaspect(toplot.ASPECT);
        if BWR
        colormap(bluewhitered)
        end
        xlabel(XNAME); ylabel(YNAME); %colorbar;
        drawnow 
        % Capture the plot as an image 
        frame = getframe(fig); 
        switch format
            case '.gif'
                im = frame2im(frame); 
                [imind,cm] = rgb2ind(im,32); 
                % Write to the GIF File 
                if in == 1 
                  imwrite(imind,cm,FILENAME,'gif', 'Loopcount',inf); 
                else 
                  imwrite(imind,cm,FILENAME,'gif','WriteMode','append', 'DelayTime',DELAY);
                end 
            case '.avi'
                writeVideo(vidfile,frame); 
            otherwise
                disp('Unknown format');
                break
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
    switch format
        case '.gif'
            disp(['Gif saved @ : ',FILENAME])
        case '.avi'
            disp(['Video saved @ : ',FILENAME])
            close(vidfile);
    end
end
end
