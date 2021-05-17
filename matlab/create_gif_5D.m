title1 = GIFNAME;
FIGDIR = BASIC.RESDIR;
GIFNAME = [FIGDIR, GIFNAME,'.gif'];

sz = size(FIELD);

% Setup figure frame
fig  = figure('Color','white','Position', [100, 100, sz(2)*400, 400]);
    for ij_ = 1:sz(2)
    subplot(100+sz(2)*10+ij_)
%         pclr = imagesc(X,Y,squeeze(FIELD(:,ij_,:,FRAMES(1))));
        pclr = imagesc(X,Y,squeeze(log(FIELD(:,ij_,:,FRAMES(1)))));
        xlabel('$k_r$');
        if ij_ == 1
            ylabel('$P$(max o. $k_z$)');
        else
            yticks([])
        end
        LEGEND = ['$|',FIELDNAME,'^{p',num2str(ij_-1),'}|$']; title(LEGEND);
    end
%     colormap gray
    axis tight manual % this ensures that getframe() returns a consistent size
    in      = 1;
    nbytes = fprintf(2,'frame %d/%d',in,numel(FIELD(1,1,1,:)));
    for n = FRAMES % loop over selected frames
%         scale = max(max(max(abs(FIELD(:,:,:,n)))));
        for ij_ = 1:sz(2)
            subplot(100+sz(2)*10+ij_)
            scale = max(max(max(abs(FIELD(:,ij_,:,n)))));
            pclr = imagesc(X,Y,squeeze(FIELD(:,ij_,:,n))/scale); caxis([0,1]);
%             pclr = imagesc(X,Y,squeeze(log(FIELD(:,ij_,:,n))));
            xlabel(XNAME);
            if ij_ == 1
                ylabel(YNAME);
            else
                yticks([])
            end
            LEGEND = ['$|',FIELDNAME,'^{p',num2str(ij_-1),'}|$']; 
            title([LEGEND,', amp = ',sprintf('%.1e',scale)]);
        end
        suptitle(['$t \approx$', sprintf('%.3d',ceil(T(n)))]);
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
        nbytes = fprintf(2,'frame %d/%d',n,numel(FIELD(1,1,1,:)));
        in = in + 1;
    end
    disp(' ')
    disp(['Gif saved @ : ',GIFNAME])

