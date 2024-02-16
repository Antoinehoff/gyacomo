function create_film(DATA,OPTIONS,format)
%% Plot options
FPS = OPTIONS.FPS; DELAY = 1/FPS;
BWR = OPTIONS.BWR; NORMALIZED = 1;
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
%%
FILENAME  = [DATA.localdir,toplot.FILENAME,format];
XNAME     = ['$',toplot.XNAME,'$'];
YNAME     = ['$',toplot.YNAME,'$'];
FIELDNAME = ['$',toplot.FIELDNAME,'$'];
FIELD     = toplot.FIELD; X = toplot.X; Y = toplot.Y;
TIME      = toplot.TIME;
switch format
    case '.avi'
        vidfile = VideoWriter(FILENAME,'Uncompressed AVI');
        vidfile.FrameRate = FPS;
        open(vidfile);  
end
% Set colormap boundaries
hmax = max(max(max(FIELD(:,:,:))));
hmin = min(min(min(FIELD(:,:,:))));
scale = -1;
flag = 0;
if hmax == hmin 
    disp('Warning : h = hmin = hmax = const')
else
% Setup figure frame
    if ~strcmp(OPTIONS.PLAN,'sx')
        if ~strcmp(OPTIONS.PLAN,'3D')
            fig  = figure('Color','white');
            pclr = pcolor(X,Y,FIELD(:,:,1)/scale); % frame plot\
            set(pclr, 'edgecolor','none'); %pbaspect(toplot.ASPECT);
        else
            fig  = figure('Color','white','Position', 1024*[0.5 0.5 0.5 1]);
            s = surface(toplot.X,toplot.Y,OPTIONS.XYZ(3)+0*toplot.X,FIELD(:,:,1)/scale); hold on;
            s.EdgeColor = 'none';
            s = surface(toplot_xz.X,OPTIONS.XYZ(2)+0*toplot_xz.X,toplot_xz.Y,toplot_xz.FIELD(:,:,1)./scale);
            s.EdgeColor = 'none';
            s = surface(OPTIONS.XYZ(1)+0*toplot_yz.X,toplot_yz.X,toplot_yz.Y,toplot_yz.FIELD(:,:,1)./scale);
            s.EdgeColor = 'none';
            zlabel('z');
            view([1 -1 0.25])          
        end        
        if BWR
            colormap(bluewhitered)
        else
            colormap(gray)
        end
        if OPTIONS.CLIMAUTO
            clim('auto')
        end
        axis tight manual % this ensures that getframe() returns a consistent size
        if toplot.INTERP
            shading interp;
        end
    else
      fig  = figure('Color','white');
      contour(toplot.X,toplot.Y,FIELD(:,:,1),128);        
    end
    in      = 1;
    nbytes = fprintf(2,'frame %d/%d',in,numel(FIELD(1,1,:)));
    for n = 1:numel(TIME) % loop over selected frames
        scale = max(max(abs(FIELD(:,:,n)))); % Scaling to normalize
        if ~strcmp(OPTIONS.PLAN,'sx')
            if ~strcmp(OPTIONS.PLAN,'3D')
                pclr = pcolor(X,Y,FIELD(:,:,n)/scale); % frame plot\
                set(pclr, 'edgecolor','none'); %pbaspect(toplot.ASPECT);
            else
                s = surface(toplot.X,toplot.Y,OPTIONS.XYZ(3)+0*toplot.X,FIELD(:,:,n)/scale); hold on;
                s.EdgeColor = 'none';
                s = surface(toplot_xz.X,OPTIONS.XYZ(2)+0*toplot_xz.X,toplot_xz.Y,toplot_xz.FIELD(:,:,n)./scale);
                s.EdgeColor = 'none';
                s = surface(OPTIONS.XYZ(1)+0*toplot_yz.X,toplot_yz.X,toplot_yz.Y,toplot_yz.FIELD(:,:,n)./scale);
                s.EdgeColor = 'none';
                zlabel('z');
                view([1 -1 0.25])          
            end
            clim([-1,1]);
            if toplot.INTERP
                shading interp; 
            end
            if BWR
                colormap(bluewhitered)
            else
                colormap(gray)
            end
            if OPTIONS.CLIMAUTO
            clim('auto')
            end
        else % show velocity distr.
           contour(toplot.X,toplot.Y,FIELD(:,:,n)/scale,128);
        end
        if ~OPTIONS.RMAXIS
            title([FIELDNAME,', $t \approx$', sprintf('%.3d',ceil(TIME(n)))...
                  ,', scaling = ',sprintf('%.1e',scale)]);
            xlabel(XNAME); ylabel(YNAME); %colorbar;
        else
            axis off
            axis tight
        end
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
