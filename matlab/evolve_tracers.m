% Options
SHOW_FILM = 0;
U_TIME   = 50; % >0 for frozen velocity at a given time, -1 for evolving field
Evolve_U = 1; % 0 for frozen velocity at a given time, 1 for evolving field
Tfin   = 50;
dt_    = 0.05;
Nstep  = ceil(Tfin/dt_);
% Init tracers
Np      = 100; %number of tracers
% color = tcolors;
color = jet(Np);
tcolors = distinguishable_colors(Np); %Their colors

Na = 10000; %length of trace

Traj_x = zeros(Np,Nstep);
Traj_y = zeros(Np,Nstep);
Disp_x = zeros(Np,Nstep);
Disp_y = zeros(Np,Nstep);

xmax = max(data.x); xmin = min(data.x);
ymax = max(data.y); ymin = min(data.y);

INIT = 'round';
switch INIT
    case 'lin'
        % Evenly distributed initial positions
        dp_ = (xmax-xmin)/(Np-1);
        Xp = linspace(xmin+dp_/2,xmax-dp_/2,Np);
        Yp = zeros(1,Np);
        Zp = zeros(Np);
    case 'round'
        % All particles arround a same point
        xc = 0; yc = 0;
        theta = rand(1,Np)*2*pi; r = 0.1;
        Xp = xc + r*cos(theta);
        Yp = yc + r*sin(theta);
        Zp = zeros(Np);
    case 'gauss'
        % normal distribution arround a point
        xc = 0; yc = 0; sgm = 1.0;
        dx = normrnd(0,sgm,[1,Np]); dx = dx - mean(dx);
        Xp = xc + dx;
        dy = normrnd(0,sgm,[1,Np]); dy = dy - mean(dy);
        Yp = yc + dy;
        Zp = zeros(Np);        
end

% position grid and velocity field
[YY_, XX_ ,ZZ_] = meshgrid(data.y,data.x,data.z);
[KX,KY] = meshgrid(data.kx,data.ky);
Ux = zeros(size(XX_));
Uy = zeros(size(XX_));
Uz = zeros(size(XX_));
ni = zeros(size(XX_));

[~,itu_] = min(abs(U_TIME-data.Ts3D));
% computing the frozen velocity field
for iz = 1:data.Nz
%     Ux(:,:,iz) = real(ifft2( 1i*KY.*(data.PHI(:,:,iz,itu_)),data.Nx,data.Ny));
%     Uy(:,:,iz) = real(ifft2(-1i*KX.*(data.PHI(:,:,iz,itu_)),data.Nx,data.Ny));
%     ni(:,:,iz) = real(ifft2(data.DENS_I(:,:,iz,itu_),data.Nx,data.Ny));
    Ux(:,:,iz) = ifourier_GENE( 1i*KY.*(data.PHI(:,:,iz,itu_)))'*sqrt(data.scale);
    Uy(:,:,iz) = ifourier_GENE(-1i*KX.*(data.PHI(:,:,iz,itu_)))'*sqrt(data.scale);
    ni(:,:,iz) = ifourier_GENE(data.DENS_I(:,:,iz,itu_))'*sqrt(data.scale);
end

%% FILM options
FPS = 30; DELAY = 1/FPS;
FORMAT = '.gif';
if SHOW_FILM
    FILENAME  = [data.localdir,'tracer_evolution',FORMAT];
    switch FORMAT
        case '.avi'
            vidfile = VideoWriter(FILENAME,'Uncompressed AVI');
            vidfile.FrameRate = FPS;
            open(vidfile);  
    end 
    fig = figure;
end

%
%Time loop
t_ = 0;
it = 1;
itu_old = 0;
nbytes = fprintf(2,'frame %d/%d',it,Nstep);
while(t_<Tfin && it <= Nstep)
   if Evolve_U
    [~,itu_] = min(abs(U_TIME+t_-data.Ts3D));
   end
    if Evolve_U && (itu_old ~= itu_)
        % updating the velocity field and density field
        for iz = 1:data.Nz
            Ux(:,:,iz) = ifourier_GENE( 1i*KY.*(data.PHI(:,:,iz,itu_)))'*sqrt(data.scale);
            Uy(:,:,iz) = ifourier_GENE(-1i*KX.*(data.PHI(:,:,iz,itu_)))'*sqrt(data.scale);
            ni(:,:,iz) = ifourier_GENE(data.DENS_I(:,:,iz,itu_))'*sqrt(data.scale);
        end
    end
    % evolve each tracer
    for ip = 1:Np
        % locate the tracer
            % find corners of the cell
            x_ = Xp(ip);
            [e_x,ixC] = min(abs(x_-data.x)); 
            if e_x == 0 % on the face
                ix0 = ixC;
                ix1 = ixC;
            elseif x_ > data.x(ixC) % right from grid point
                ix0 = ixC;
                ix1 = ixC+1;
            else % left
                ix0 = ixC-1;  
                ix1 = ixC; 
            end
            y_ = Yp(ip);
            [e_y,iyC] = min(abs(y_-data.y));
            if e_y == 0 % on the face
                iy0 = iyC;
                iy1 = iyC;
            elseif y_ > data.y(iyC) % above
                iy0 = iyC;
                iy1 = iyC+1;
            else % under
                iy0 = iyC-1;  
                iy1 = iyC; 
            end
            z_ = Zp(ip,1);
            [e_z,izC] = min(abs(z_-data.z));
            if e_z == 0 % on the face
                iz0 = izC;
                iz1 = izC;
            elseif z_ > data.z(izC) % before
                iz0 = izC;
                iz1 = izC+1;
            else % behind
                iz0 = izC-1;  
                iz1 = izC; 
            end
            x0   = data.x(ix0); x1 = data.x(ix1); %left right
            y0   = data.y(iy0); y1 = data.y(iy1); %down top
            z0   = data.z(iz0); z1 = data.z(iz1); %back front
            if(e_x > 0)
                ai__ = (x_ - x0)/(x1-x0); % interp coeff x
            else
                ai__ = 0;
            end
            if(e_y > 0)
                a_i_ = (y_ - y0)/(y0-y1); % interp coeff y
            else
                a_i_ = 0;
            end
            if(e_z > 0)
                a__i = (z_ - z0)/(z1-z0); % interp coeff z
            else
                a__i = 0;
            end
        % interp velocity and density
            %velocity values
            u000 = [Ux(ix0,iy0,iz0) Uy(ix0,iy0,iz0) ni(ix0,iy0,iz0)];
            u001 = [Ux(ix0,iy0,iz1) Uy(ix0,iy0,iz1) ni(ix0,iy0,iz1)];
            u010 = [Ux(ix0,iy1,iz0) Uy(ix0,iy1,iz0) ni(ix0,iy1,iz0)];
            u011 = [Ux(ix0,iy1,iz1) Uy(ix0,iy1,iz1) ni(ix0,iy1,iz1)];
            u100 = [Ux(ix1,iy0,iz0) Uy(ix1,iy0,iz0) ni(ix1,iy0,iz0)];
            u101 = [Ux(ix1,iy0,iz1) Uy(ix1,iy0,iz1) ni(ix1,iy0,iz1)];
            u110 = [Ux(ix1,iy1,iz0) Uy(ix1,iy1,iz0) ni(ix1,iy1,iz0)];
            u111 = [Ux(ix1,iy1,iz1) Uy(ix1,iy1,iz1) ni(ix1,iy1,iz1)];
            %linear interpolation
            linterp = @(x1,x2,a) (1-a)*x1 + a*x2;
            
            u_00  =  linterp(u000,u100,ai__);
            u_10  =  linterp(u010,u110,ai__);
            u__0  =  linterp(u_00,u_10,a_i_);
            
            u_01  =  linterp(u001,u101,ai__);
            u_11  =  linterp(u011,u111,ai__);
            u__1  =  linterp(u_01,u_11,a_i_);
            
            u___  =  linterp(u__0,u__1,a__i);

%             push the particle
%             q = sign(-u___(3));
            q = -u___(3);
%             q =1;
            x_ = x_ + dt_*u___(1)*q;
            y_ = y_ + dt_*u___(2)*q;
                
        % apply periodic boundary conditions
        if(x_ > xmax)
            x_ = xmin + (x_ - xmax);
        elseif(x_ < xmin)
            x_ = xmax + (x_ - xmin);
        end
        if(y_ > ymax)
            y_ = ymin + (y_ - ymax);
        elseif(y_ < ymin)
            y_ = ymax + (y_ - ymin);
        end
        % store data
            % Trajectory
            Traj_x(ip,it) = x_;
            Traj_y(ip,it) = y_;  
            % Displacement
            Disp_x(ip,it) = Disp_x(ip,it) + dt_*u___(1)*q;
            Disp_y(ip,it) = Disp_y(ip,it) + dt_*u___(2)*q;        
        % update position
        Xp(ip) = x_; Yp(ip) = y_;
    end
    %% Movie
    if SHOW_FILM && (~Evolve_U || (itu_old ~= itu_))
    % updating the velocity field
        clf(fig);
        F2P = ifourier_GENE(data.PHI(:,:,iz,itu_))';
        scale = max(max(abs(F2P))); % Scaling to normalize
        pclr = pcolor(XX_,YY_,F2P/scale); 
        colormap(bluewhitered);
        set(pclr, 'edgecolor','none'); hold on; caxis([-2,2]); shading interp
        for ip = 1:Np
            ia0 = max(1,it-Na);
            plot(Traj_x(ip,ia0:it),Traj_y(ip,ia0:it),'.','Color',color(ip,:)); hold on
        end
        for ip = 1:Np
            plot(Traj_x(ip, 1),Traj_y(ip, 1),'xk'); hold on
            plot(Traj_x(ip,it),Traj_y(ip,it),'ok','MarkerFaceColor',color(ip,:)); hold on
        end
        title(['$t \approx$', sprintf('%.3d',ceil(data.Ts3D(itu_)))]);
        axis equal
        xlim([xmin xmax]); ylim([ymin ymax]);
        drawnow
        % Capture the plot as an image 
        frame = getframe(fig); 
        switch FORMAT
            case '.gif'
                im = frame2im(frame); 
                [imind,cm] = rgb2ind(im,32); 
                % Write to the GIF File 
                if it == 1 
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
    end
    t_ = t_ + dt_; it = it + 1; itu_old = itu_;
    % terminal info
    while nbytes > 0
      fprintf('\b')
      nbytes = nbytes - 1;
    end
    nbytes = fprintf(2,'frame %d/%d',it,Nstep);
end
disp(' ')
switch FORMAT
    case '.gif'
        disp(['Gif saved @ : ',FILENAME])
    case '.avi'
        disp(['Video saved @ : ',FILENAME])
        close(vidfile);
end

Nt = it;
%% Plot trajectories and statistics
xtot = Disp_x;
ytot = Disp_y;
% aggregate the displacements
for iq = 1:Np
    x_ = 0; y_ = 0;
    for iu = 1:Nt-1
            x_          = x_ + Disp_x(iq,iu);
            y_          = y_ + Disp_y(iq,iu);
            xtot(iq,iu) = x_;
            ytot(iq,iu) = y_;
    end
end
ma = @(x) x;%movmean(x,100);
time_ = linspace(U_TIME,U_TIME+Tfin,it-1);
figure;
subplot(221);
for ip = 1:Np
%     plot(ma(time_),ma(Traj_x(ip,:)),'-','Color',tcolors(ip,:)); hold on
    plot(ma(time_),xtot(ip,:),'-','Color',tcolors(ip,:)); hold on
%     plot(ma(time_),ma(Disp_x(ip,:)),'-','Color',tcolors(ip,:)); hold on
end
ylabel('$x_p$');
xlim(U_TIME + [0 Tfin]);

subplot(222);
    itf = floor(Nt/2); %fit end time
    % x^2 displacement
    plot(time_,mean(xtot.^2,1),'DisplayName','$\langle x.^2\rangle_p$'); hold on
    fit = polyfit(time_(1:itf),mean(xtot(:,1:itf).^2,1),1);
    plot(time_,fit(1)*time_+fit(2),'--k'); hold on
    ylabel('$\langle x^2 \rangle_p$');

%     % y^2 displacement
%     fit = polyfit(time_(1:itf),mean(ytot(:,1:itf).^2,1),1);
%     plot(time_,fit(1)*time_+fit(2),'--k','DisplayName',['$\alpha=',num2str(fit(1)),'$']); hold on
%     plot(time_,mean(ytot.^2,1),'DisplayName','$\langle y.^2\rangle_p$'); 
%     
%     % r^2 displacement
%     fit = polyfit(time_(1:itf),mean(xtot(:,1:itf).^2+ytot(:,1:itf).^2,1),1);
%     plot(time_,fit(1)*time_+fit(2),'--k','DisplayName',['$\alpha=',num2str(fit(1)),'$']); hold on
%     plot(time_,mean(xtot.^2+ytot.^2,1),'DisplayName','$\langle r.^2\rangle_p$'); 
%     ylabel('$\langle x^2 \rangle_p$');
%     xlim(U_TIME + [0 Tfin]);

subplot(223);
for ip = 1:Np
%     plot(ma(time_),ma(Traj_y(ip,:)),'-','Color',tcolors(ip,:)); hold on
    plot(ma(time_),ytot(ip,:),'-','Color',tcolors(ip,:)); hold on
%     plot(ma(time_),ma(Disp_y(ip,:)),'-','Color',tcolors(ip,:)); hold on
end
xlabel('time');
ylabel('$y_p$');
xlim(U_TIME + [0 Tfin]);

subplot(224);
histogram(xtot(:,1),20); hold on
histogram(xtot(:,end),20)
xlabel('position');
ylabel('$n$');

