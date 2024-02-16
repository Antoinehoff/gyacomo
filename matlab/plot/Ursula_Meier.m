function [filename] = Ursula_Meier(MOVIE)
    FPS = 16; BWR = 0; ncolor_level = 128; format = '.gif'; SAT = 0.75;
    INTERP = 1;
    if isfield(MOVIE,'FPS') 
        FPS = MOVIE.FPS;
    end
    if isfield(MOVIE,'BWR') 
        BWR = MOVIE.BWR;
    end
    if isfield(MOVIE,'format') 
        format = MOVIE.format;
    end
    if isfield(MOVIE,'SAT') 
        SAT = MOVIE.SAT;
    end
    if isfield(MOVIE,'INTERP') 
        INTERP = MOVIE.INTERP;
    end
    DELAY = 1/MOVIE.FPS;
    FILENAME = MOVIE.FILENAME; FIELDNAME = MOVIE.FIELDNAME;
    X = MOVIE.X; Y = MOVIE.Y; T = MOVIE.T; F = MOVIE.F;
    XNAME = MOVIE.XNAME; YNAME = MOVIE.YNAME;
    switch format
        case '.avi'
            vidfile = VideoWriter(FILENAME,'Uncompressed AVI');
            vidfile.FrameRate = FPS;
            open(vidfile);  
    end
    % Set colormap boundaries
    hmax = max(max(max(F(:,:,:))));
    hmin = min(min(min(F(:,:,:))));
    if hmax == hmin 
        disp('Warning : h = hmin = hmax = const')
    else
    % Setup figure frame
    fig  = figure('Color','white');%,'Position', toplot.DIMENSIONS.*[0.5 0.5 1.0 1]);
        pcolor(X,Y,F(:,:,1)); % to set up
        if BWR
            colormap(bluewhitered)
        else
            colormap(gray)
        end
        clim(SAT*[-1 1])
        axis equal
        axis tight % this ensures that getframe() returns a consistent size
        if INTERP
            shading interp;
        end
        in      = 1;
        nbytes = fprintf(2,'frame %d/%d',in,numel(F(1,1,:)));
        for n = 1:numel(T) % loop over selected frames
            scale = max(max(abs(F(:,:,n)))); % Scaling to normalize
            pclr = pcolor(X,Y,F(:,:,n)/scale); % frame plot\
            clim([-1,1]);
            if INTERP
                shading interp; 
            end
            set(pclr, 'edgecolor','none'); %pbaspect(toplot.ASPECT);
            if BWR
                colormap(bluewhitered)
            else
                colormap(gray)
            end
            clim(SAT*[-1 1])
            title([FIELDNAME,', $t \approx$', sprintf('%.3d',ceil(T(n)))...
                  ,', scaling = ',sprintf('%.1e',scale)]);
    
            xlabel(XNAME); ylabel(YNAME); %colorbar;
            axis tight equal
            drawnow 
            % Capture the plot as an image 
            frame = getframe(fig); 
            switch format
                case '.gif'
                    im = frame2im(frame); 
                    [imind,cm] = rgb2ind(im,ncolor_level); 
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
            nbytes = fprintf(2,'frame %d/%d',n,numel(F(1,1,:)));
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