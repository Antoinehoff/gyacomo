function [ ] = create_gif(x, y, t, h, BASIC, GRID, MODEL, delay, GIFNAME, saturate)
title1 = GIFNAME;
FIGDIR = ['../results/', BASIC.SIMID,'/'];
if ~exist(FIGDIR, 'dir')
   mkdir(FIGDIR)
end

if     MODEL.CO == -1; CONAME = 'FC';
elseif MODEL.CO == -2; CONAME = 'DC';
elseif MODEL.CO ==  0; CONAME = 'LB'; end;

GIFNAME = [GIFNAME,'_Pe_',num2str(GRID.pmaxe),'_Je_',num2str(GRID.jmaxe),...
    '_Pi_',num2str(GRID.pmaxi),'_Ji_',num2str(GRID.jmaxi),...
    '_etan_',num2str(MODEL.eta_n),'_etaB_',num2str(MODEL.eta_B),'_etaT_',...
    num2str(MODEL.eta_T),'_nu_',num2str(MODEL.nu),'_',CONAME];
GIFNAME = [FIGDIR, GIFNAME,'.gif'];

% Set colormap boundaries
hmax = max(max(max(h)));
hmin = min(min(min(h)));

flag = 0;
if hmax == hmin 
    disp('Warning : h = hmin = hmax = const')
else
    % Setup figure frame
    fig  = figure;
    pcolor(x,y,h(:,:,1)); % to set up
    colormap jet
    colorbar
    if not(flag) && saturate
        caxis([hmin hmax])
    end
    axis tight manual % this ensures that getframe() returns a consistent size
    hold on
    
    n      = 0;
    nbytes = fprintf(2,'frame %d/%d',n,numel(h(1,1,:)));
    for n = 1:numel(h(1,1,:)) % time loop
        pclr = pcolor(x,y,h(:,:,n)); % frame plot
        set(pclr, 'edgecolor','none')
        title([title1,', $t \approx$', sprintf('%.3d',ceil(t(n)))])
        drawnow 
        % Capture the plot as an image 
        frame = getframe(fig); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,64); 
        % Write to the GIF File 
        if n == 1 
          imwrite(imind,cm,GIFNAME,'gif', 'Loopcount',inf); 
        else 
          imwrite(imind,cm,GIFNAME,'gif','WriteMode','append', 'DelayTime',delay); 
        end 
        % terminal info
        while nbytes > 0
          fprintf('\b')
          nbytes = nbytes - 1;
        end
        nbytes = fprintf(2,'frame %d/%d',n,numel(h(1,1,:)));
    end

    disp(['Gif saved @ : ',GIFNAME])
end
end

