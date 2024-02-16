%% This script is made for illustrating different magnetic geometry
% it is not bugfree but work pretty well for illustration purpose of the
% magnetic flux surface of a tokamak and for the particle trajectory in it.
% ! The particle trajectory is not correct (does not work for Z-pinch for
% ! example) but is good enough for illustation.
gyacomodir = pwd; gyacomodir = gyacomodir(1:end-2); % get code directory
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add
default_plots_options
% Torus flux tube geometry

Nx     = 128;
Ny     = 128;
Nz     = 64;
Nturns = 1.0;
lr     = 0.05; %larmor radius
Wc     = 750;  %~cyclotron freq
vd     = 0.25; %~drift velocty
Gx      = linspace(-60,60,Nx)*0.001;
y      = linspace(-60,60,Ny)*0.001;
FIELD  = ones(Nx,Ny,Nz);
z      = linspace(-Nturns*pi,Nturns*pi,Nz);
N_field_lines    = 10;
N_magn_flux_surf = 1;
openangle        = pi/3;
% phi0 = openangle/2; phi1= 2*pi-openangle/2;
phi0 = 0; phi1= 2*pi-1.3;
PLT_FTUBE = 0;
PLT_BASIS = 0;
PLT_DATA  = 0;
R0     = 1; % Torus major radius
Z0     = 0;
xpoint = 0;
select = 1;
switch select
    case 1
    % Z-pinch
    rho    = 0.50; drho = 0.1;% Torus minor radius
    eps    = 0.3;
    q0     = 0.000001; % Flux tube safety factor
    shear  = 0.8;
    kappa  = 1.0;
    s_kappa= 0.0;
    delta  = 0.0;
    s_delta= 0.0;
    zeta   = 0.0;
    s_zeta = 0.0;
    case 2
    % CBC
    rho    = 0.50; drho = 0.1;% Torus minor radius
    eps    = 0.18;
    q0     = 1.4; % Flux tube safety factor
    shear  = 0.8;
    kappa  = 1.0;
    s_kappa= 0.0;
    delta  = 0.0;
    s_delta= 0.0;
    zeta   = 0.0;
    s_zeta = 0.0;
    case 3
    % DIII-D edge
    rho    = 0.90; drho = 0.1;% Torus minor radius
    eps    = 0.3;
    q0     = 4.8; % Flux tube safety factor
    shear  = 2.5;
    kappa  = 1.55;
    s_kappa= 1.0;
    delta  = 0.4;
    s_delta= 1.5;
    zeta   = 0.1;
    s_zeta = 0.0;
    case 4 
    % Poloidal DIII-D
    N_magn_flux_surf = 5;
    phi0 = 0; phi1= pi/48;
    N_field_lines    = 0;
    rho    = 0.50; drho = 0.1;% Torus minor radius
    eps    = 0.3;
    q0     = 4.8; % Flux tube safety factor
    shear  = 2.5;
    kappa  = 1.55;
    s_kappa= 2.0;
    delta  = 0.3;
    s_delta= 3.0;
    zeta   = 0.1;
    s_zeta = 0.0;
    case 5
    % play with DIII-D edge
    rho    = 0.90; drho = 0.1;% Torus minor radius
    eps    = 0.3;
    q0     = 1.8; % Flux tube safety factor
    % Nturns = max(1,q0);
    shear  = 2.5;
    kappa  = 1.5;
    s_kappa= 1.0;
    delta  = 0.3;
    s_delta= 0.0;
    zeta   = 0.0;
    s_zeta = 0.0;
end
% Nturns = 1/q0;
Nptor  = 64;
% Toroidal angle for the flux surf
% phi    = linspace(0+openangle/2, 2*pi-openangle/2, Nptor);
phi    = linspace(phi0,phi1, Nptor);
theta  = linspace(-pi, pi, Nptor)   ; % Poloidal angle
[p, t] = meshgrid(phi, theta);
Tgeom  = 1;
DIMENSIONS = [600 600 1200 600];
rho = rho*eps; drho = drho*eps;
iota = 1/q0;
% field line coordinates
Xfl = @(z) (R0+rho*cos(z+asin(delta*sin(z)))).*cos(q0*z);
Yfl = @(z) (R0+rho*cos(z+asin(delta*sin(z)))).*sin(q0*z);
Zfl = @(z)  Z0+kappa*rho*sin(z+zeta*sin(2*z)+xpoint*cos(2*z));
Rvec= @(z) [Xfl(z); Yfl(z); Zfl(z)];
% xvec shearless
% xX  = @(z) (Xfl(z)-R0*cos(q0*z))./sqrt((Xfl(z)-R0*cos(q0*z)).^2+(Yfl(z)-R0*sin(q0*z)).^2+Zfl(z).^2);
% xY  = @(z) (Yfl(z)-R0*sin(q0*z))./sqrt((Xfl(z)-R0*cos(q0*z)).^2+(Yfl(z)-R0*sin(q0*z)).^2+Zfl(z).^2);
% xZ  = @(z)              Zfl(z)./sqrt((Xfl(z)-R0*cos(q0*z)).^2+(Yfl(z)-R0*sin(q0*z)).^2+Zfl(z).^2);
xX = @(z) cos(z+asin(delta*sin(z))).*cos(q0*z)./sqrt(cos(z+asin(delta*sin(z))).^2+kappa^2*sin(z+zeta*sin(2*z)).^2);
xY = @(z) cos(z+asin(delta*sin(z))).*sin(q0*z)./sqrt(cos(z+asin(delta*sin(z))).^2+kappa^2*sin(z+zeta*sin(2*z)).^2);
xZ = @(z) kappa*sin(z+zeta*sin(2*z))./sqrt(cos(z+asin(delta*sin(z))).^2+kappa^2*sin(z+zeta*sin(2*z)).^2);
xvec= @(z) [xX(z); xY(z); xZ(z)];
% bvec shearless
bX  = @(z) Tgeom*(rho*cos(z).*cos(q0*z) - q0*Yfl(z))./sqrt(Tgeom*(rho*cos(z).*cos(q0*z) - q0*Yfl(z)).^2+(rho*cos(z).*sin(q0*z) + q0*Xfl(z)).^2+(rho*cos(z)).^2);
bY  = @(z)       (rho*cos(z).*sin(q0*z) + q0*Xfl(z))./sqrt(Tgeom*(rho*cos(z).*cos(q0*z) - q0*Yfl(z)).^2+(rho*cos(z).*sin(q0*z) + q0*Xfl(z)).^2+(rho*cos(z)).^2);
bZ  = @(z)                              rho*cos(z)./sqrt(Tgeom*(rho*cos(z).*cos(q0*z) - q0*Yfl(z)).^2+(rho*cos(z).*sin(q0*z) + q0*Xfl(z)).^2+(rho*cos(z)).^2);
bvec= @(z) [bX(z); bY(z); bZ(z)];
% yvec = b times x
yX  = @(z) bY(z).*xZ(z)-bZ(z).*xY(z)./sqrt((bY(z).*xZ(z)-bZ(z).*xY(z)).^2+(bZ(z).*xX(z)-bX(z).*xZ(z)).^2+(bX(z).*xY(z)-bY(z).*xX(z)).^2);
yY  = @(z) bZ(z).*xX(z)-bX(z).*xZ(z)./sqrt((bY(z).*xZ(z)-bZ(z).*xY(z)).^2+(bZ(z).*xX(z)-bX(z).*xZ(z)).^2+(bX(z).*xY(z)-bY(z).*xX(z)).^2);
yZ  = @(z) bX(z).*xY(z)-bY(z).*xX(z)./sqrt((bY(z).*xZ(z)-bZ(z).*xY(z)).^2+(bZ(z).*xX(z)-bX(z).*xZ(z)).^2+(bX(z).*xY(z)-bY(z).*xX(z)).^2);
yvec= @(z) [yX(z); yY(z); yZ(z)];

scale = 0.10;

figure; set(gcf, 'Position',  DIMENSIONS)

    %plot magnetic geometry
    if N_magn_flux_surf > 0
        dr_array = drho*(-(N_magn_flux_surf-1):(N_magn_flux_surf-1));
        Nsurf    = numel(dr_array);
        clrs     = cool(Nsurf);
        % clrs     = clrs(end:-1:1,:);
        for is = 1:Nsurf
            [xt, yt, zt] = ...
                mag_flux_surf...
                    (t,p,R0,Z0,rho,dr_array(is),kappa,s_kappa,delta,s_delta,zeta,s_zeta,xpoint);
            magnetic_topo=surf(xt, yt, zt); hold on;
            alpha = 0.7*(1-(is/Nsurf)^2)+0.5*(Nsurf==1);
            set(magnetic_topo,...
                'edgecolor','none',...
                'facecolor',clrs(is,:),...
                'FaceAlpha',alpha);
        end
        H = light;
    end
    %plot field lines
    if N_field_lines > 0
        % theta  = linspace(-Nturns*pi, Nturns*pi, 1024)   ; % Poloidal angle
        theta  = linspace(0, Nturns*2*pi, 1024)   ; % Poloidal angle
        vecfl = Rvec(theta);
        plot3(vecfl(1,:),vecfl(2,:),vecfl(3,:),'-b'); hold on;
        % Multiple field lines
        x_ = []; y_ = []; z_ = [];
        for ifl = 1:N_field_lines-1
            % rotation for multiple field lines
            t  = q0*pi*ifl;
            % x_ = [x_ cos(t)*vecfl(1,:) - sin(t)*vecfl(2,:)];
            % y_ = [y_ sin(t)*vecfl(1,:) + cos(t)*vecfl(2,:)];
            % z_ = [z_ vecfl(3,:)];
            x_ = cos(t)*vecfl(1,:) - sin(t)*vecfl(2,:);
            y_ = sin(t)*vecfl(1,:) + cos(t)*vecfl(2,:);
            z_ = vecfl(3,:);
            % plot3(x_,y_,z_,'-','color',[1.0 0.6 0.6]*0.8); hold on;
            plot3(x_,y_,z_,'-','color',[0.0 0.5 1.0]*0.8,'LineWidth',1); hold on;
        end
        %plot3(x_,y_,z_,'-','color',[1.0 0.6 0.6]*0.8); hold on;
    end
    %plot fluxe tube
    if PLT_FTUBE
        theta  = linspace(-Nturns*pi, Nturns*pi, 64)    ; % Poloidal angle
        %store the shifts in an order (top left to bottom right)
        s_x = rhostar*[Gx(1) Gx(end) Gx(1)   Gx(end)]; 
        s_y = rhostar*[y(1) y(1)   y(end) y(end)];
        for i_ = 1:4
        vx_ = Xfl(theta) + s_x(i_)*xX(theta) + s_y(i_)*yX(theta);
        vy_ = Yfl(theta) + s_x(i_)*xY(theta) + s_y(i_)*yY(theta);
        vz_ = Zfl(theta) + s_x(i_)*xZ(theta) + s_y(i_)*yZ(theta);
        plot3(vx_,vy_,vz_,'-','color',[1.0 0.6 0.6]*0.8,'linewidth',1.5); hold on;
        end
    end
    %plot vector basis
    if PLT_BASIS
        theta  = z   ; % Poloidal angle
        plot3(Xfl(theta),Yfl(theta),Zfl(theta),'ok'); hold on;
        quiver3(Xfl(theta),Yfl(theta),Zfl(theta),scale*xX(theta),scale*xY(theta),scale*xZ(theta),0,'r');
        quiver3(Xfl(theta),Yfl(theta),Zfl(theta),scale*yX(theta),scale*yY(theta),scale*yZ(theta),0,'g');
        quiver3(Xfl(theta),Yfl(theta),Zfl(theta),scale*bX(theta),scale*bY(theta),scale*bZ(theta),0,'b');
    end
    % Particle trajectory
    if 1
        % Particle gyration
        Gx = @(z) lr*cos(z).*cos(Wc*z);
        Gy = @(z) lr*sin(z).*cos(Wc*z);
        Gz = @(z) lr*sin(Wc*z);
        % This is the "real" helicoidal trajectory
        theta  = linspace(0, Nturns*2*pi, 10000)   ; % Poloidal angle
        dtheta = theta(2)-theta(1);
        dX = 0*theta; dY = dX; dZ = dX; % init drifts
        for i = 2:numel(theta)
            th_old = theta(i-1);
            th_    = theta(i);
            % find magnetic displacement vector
            bx_= Xfl(th_)-Xfl(th_old);
            by_= Yfl(th_)-Yfl(th_old);
            bz_= Zfl(th_)-Zfl(th_old);
            % find binormal
            yx_= xY(th_).*bz_ - xZ(th_).*by_;
            yy_= xZ(th_).*bx_ - xX(th_).*bz_;
            yz_= xX(th_).*by_ - xY(th_).*bx_;
            % normalize
            yx_=yx_/sqrt(yx_^2+yy_^2+yz_^2);
            yy_=yy_/sqrt(yx_^2+yy_^2+yz_^2);
            yz_=yz_/sqrt(yx_^2+yy_^2+yz_^2);
            % add drift
            dX(i) = dX(i-1) - vd*yx_*dtheta;
            dY(i) = dY(i-1) - vd*yy_*dtheta;
            dZ(i) = dZ(i-1) - vd*yz_*dtheta;
        end
        % add gyration
        dX = dX +  Gx(theta);
        dY = dY +  Gy(theta);
        dZ = dZ +  Gz(theta);
        % add fieldline 
        Tx = Xfl(theta)+dX;
        Ty = Yfl(theta)+dY;
        Tz = Zfl(theta)+dZ;
        plot3(Tx,Ty,Tz,'-r','LineWidth',1.2); hold on;
        plot3(Tx(1),Ty(1),Tz(1),'xk','MarkerSize',10); hold on;
        plot3(Tx(end),Ty(end),Tz(end),'ok','MarkerSize',8,'MarkerFaceColor','k'); hold on;
    end
    xlabel('X');ylabel('Y');zlabel('Z');
    %
    axis equal
    view([1,-2,1])
    grid on
    axis off

    function [x, y, z] = mag_flux_surf(theta,phi,R0,Z0,rho,dr,kap,s_k,del,s_d,zet,s_z,xp)
rho = rho + dr;
kap = kap + s_k*dr;
del = del + s_d*dr;
zet = zet + s_z*dr;
x = (R0 + rho.*cos(theta + asin(del*sin(theta)))) .* cos(phi);
y = (R0 + rho.*cos(theta + asin(del*sin(theta)))) .* sin(phi);
% z =  Z0 + kap*rho.*sin(theta + zet*sin(2*theta) + xp*cos(2*theta));
z =  Z0 + kap*rho.*(sin(theta + zet*sin(2*theta)) + xp*sin(theta+xp*cos(theta)));
    end
