
GIFNAME = ['Nipj_kr',sprintf('_%.2d',JOBNUM)]; INTERP = 0;
plt = @(x) squeeze(max((abs(x)),[],4));
FIELD = plt(Nipj); X = kr'; Y = Pi'; T = Ts5D; FRAMES = FRAMES_5D;
FIELDNAME = 'N_i'; XNAME = '$k_{max}\rho_s$'; YNAME = '$P$';

title1 = GIFNAME;
FIGDIR = BASIC.RESDIR;
GIFNAME = [FIGDIR, GIFNAME,'.gif'];

sz = size(FIELD);

% Setup figure frame
fig  = figure('Color','white','Position', [100, 100, sz(2)*400, 400]);
    for ij_ = 1:sz(2)
    subplot(100+sz(2)*10+ij_)
        pclr = imagesc(X,Y,squeeze(FIELD(:,ij_,:,1)));
        xlabel('$k_r$');
        if ij_ == 1
            ylabel('$P$(max o. $k_z$)');
        else
            yticks([])
        end
        LEGEND = ['$|',FIELDNAME,'^{p',num2str(ij_-1),'}|$']; title(LEGEND);
    end
    colormap gray
    axis tight manual % this ensures that getframe() returns a consistent size
    in      = 1;
    nbytes = fprintf(2,'frame %d/%d',in,numel(FIELD(1,1,1,:)));
    for n = FRAMES % loop over selected frames
        scale = max(max(max(abs(FIELD(:,:,:,n)))));
        for ij_ = 1:sz(2)
            subplot(100+sz(2)*10+ij_)
            pclr = imagesc(X,Y,squeeze(FIELD(:,ij_,:,n))/scale);
            xlabel(XNAME);
            if ij_ == 1
                ylabel(YNAME);
            else
                yticks([])
            end
            LEGEND = ['$|',FIELDNAME,'^{p',num2str(ij_-1),'}|$']; title(LEGEND);
        end
        suptitle(['$t \approx$', sprintf('%.3d',ceil(T(n)))...
            ,', scaling = ',sprintf('%.1e',scale)]);
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

