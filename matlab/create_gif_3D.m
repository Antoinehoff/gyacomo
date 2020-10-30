title1 = GIFNAME;
FIGDIR = BASIC.RESDIR;
GIFNAME = [FIGDIR, GIFNAME,'.gif'];

% Setup figure frame
fig  = figure('Color','white','Position', [100, 100, 400, 400]);
    plt = @(x) squeeze(reshape(x(:,:,1),[],1));
    scale_x = max(abs(plt(X(:,:,1))));
    scale_y = max(abs(plt(Y(:,:,1))));
    scale_z = max(abs(plt(Z(:,:,1))));
    plot3(plt(X)/scale_x,plt(Y)/scale_y,plt(Z)/scale_z,'.k','MarkerSize',MARKERSIZE);
    view(VIEW);
    colormap jet
    axis tight manual % this ensures that getframe() returns a consistent size
    in      = 1;
    nbytes = fprintf(2,'frame %d/%d',in,numel(FRAMES));
    for n = FRAMES % loop over selected frames
        plt = @(x) squeeze(reshape(x(:,:,n),[],1));
        scale_x = max(abs(plt(X)));
        scale_y = max(abs(plt(Y)));
        scale_z = max(abs(plt(Z)));
        plot3(plt(X)/scale_x,plt(Y)/scale_y,plt(Z)/scale_z,'.k','MarkerSize',MARKERSIZE);
        view(VIEW);
        xlabel(XNAME); ylabel(YNAME); zlabel(ZNAME);
        title(['$t \approx$', sprintf('%.3d',ceil(T(n)))...
            ,', scaling = ',sprintf('%.1e',scale_x)]);
        grid on;
        drawnow 
        % Capture the plot as an image 
        frame = getframe(fig); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,32); 
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
        nbytes = fprintf(2,'frame %d/%d',n,numel(FRAMES));
        in = in + 1;
    end
    disp(' ')
    disp(['Gif saved @ : ',GIFNAME])

