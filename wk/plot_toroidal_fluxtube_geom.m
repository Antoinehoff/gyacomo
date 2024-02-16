gyacomodir = pwd; gyacomodir = gyacomodir(1:end-2); % get code directory
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add
default_plots_options
% Torus flux tube geometry

Nx     = 128;
Ny     = 128;
Nturns = 1;
Nz     = 128;
x      = linspace(-60,60,Nx)*0.001;
y      = linspace(-60,60,Ny)*0.001;
FIELD  = ones(Nx,Ny,Nz);
z      = linspace(-Nturns*pi,Nturns*pi,Nz);
N_field_lines    = 2;
N_magn_flux_surf = 1;
openangle        = pi/3;
phi0 = openangle/2; phi1= 2*pi-openangle/2;
PLT_FTUBE = 0;
PLT_BASIS = 0;
PLT_DATA  = 0;
R0     = 1; % Torus major radius
Z0     = 0;
xpoint = 0;
select = 3;
switch select
    case 1
    % CBC
    rho    = 0.3; drho = 0.1;% Torus minor radius
    eps    = 0.3;
    q0     = 1.4; % Flux tube safety factor
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
    % Fake x-point DIII-D
    N_magn_flux_surf = 2;
    phi0 = 0; phi1= pi/48;
    N_field_lines    = 0;
    rho    = 0.40; drho = 0.1;% Torus minor radius
    eps    = 0.3;
    q0     = 4.8; % Flux tube safety factor
    shear  = 2.5;
    kappa  = 1.0;
    s_kappa= 0.0;
    delta  = 0.0;
    s_delta= 0.0;
    zeta   = 0.0;
    s_zeta = 0.0;
    xpoint = 0.5;
    case 6 
    % play with DIII-D edge
    rho    = 0.90; drho = 0.1;% Torus minor radius
    eps    = 0.3;
    q0     = 1.0; % Flux tube safety factor
    Nturns = max(1,q0);
    shear  = 2.5;
    kappa  = 0.5;
    s_kappa= 1.0;
    delta  = 0.0;
    s_delta= 0.0;
    zeta   =-1.0;
    s_zeta = 0.0;
end

Nptor  = 64;
% Toroidal angle for the flux surf
% phi    = linspace(0+openangle/2, 2*pi-openangle/2, Nptor);
phi    = linspace(phi0,phi1, Nptor);
theta  = linspace(-pi, pi, Nptor)   ; % Poloidal angle
[p, t] = meshgrid(phi, theta);
Tgeom  = 1;

DIMENSIONS = [600 600 1200 600];
rho = rho*eps; drho = drho*eps;
% field line coordinates
Xfl = @(z) (R0+rho*cos(z+asin(delta*sin(z)))).*cos(q0*z);
Yfl = @(z) (R0+rho*cos(z+asin(delta*sin(z)))).*sin(q0*z);
Zfl = @(z)  Z0+kappa*rho*sin(z+zeta*sin(2*z)+xpoint*cos(2*z));
Rvec= @(z) [Xfl(z); Yfl(z); Zfl(z)];
% xvec shearless
xX  = @(z) (Xfl(z)-R0*cos(q0*z))./sqrt((Xfl(z)-R0*cos(q0*z)).^2+(Yfl(z)-R0*sin(q0*z)).^2+Zfl(z).^2);
xY  = @(z) (Yfl(z)-R0*sin(q0*z))./sqrt((Xfl(z)-R0*cos(q0*z)).^2+(Yfl(z)-R0*sin(q0*z)).^2+Zfl(z).^2);
xZ  = @(z)              Zfl(z)./sqrt((Xfl(z)-R0*cos(q0*z)).^2+(Yfl(z)-R0*sin(q0*z)).^2+Zfl(z).^2);
xvec= @(z) [xX(z); xY(z); xZ(z)];
% bvec
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

% Plot plane result
OPTIONS.POLARPLOT = 0;
OPTIONS.PLAN      = 'xy';
rhostar           = 0.1;
[X,Y]             = meshgrid(x,y);
max_              = 0;

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
        theta  = linspace(-Nturns*pi, Nturns*pi, Nturns*256)   ; % Poloidal angle
        vecfl = Rvec(theta);
        plot3(vecfl(1,:),vecfl(2,:),vecfl(3,:),'-b'); hold on;
        % Multiple field lines
        x_ = []; y_ = []; z_ = [];
        for ifl = 1:N_field_lines-1
            % rotation for multiple field lines
            t  = q0*pi*ifl;
            x_ = [x_ cos(t)*vecfl(1,:) - sin(t)*vecfl(2,:)];
            y_ = [y_ sin(t)*vecfl(1,:) + cos(t)*vecfl(2,:)];
            z_ = [z_ vecfl(3,:)];
        end
        plot3(x_,y_,z_,'-','color',[1.0 0.6 0.6]*0.8); hold on;
    end
    %plot fluxe tube
    if PLT_FTUBE
        theta  = linspace(-Nturns*pi, Nturns*pi, 64)    ; % Poloidal angle
        %store the shifts in an order (top left to bottom right)
        s_x = rhostar*[x(1) x(end) x(1)   x(end)]; 
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
    xlabel('X');ylabel('Y');zlabel('Z');
    %Plot time dependent data
    if PLT_DATA
        for iz = 1:Nz
            z_ = z(iz);    
            Xp = Xfl(z_) + X*xX(z_) + Y*yX(z_);
            Yp = Yfl(z_) + X*xY(z_) + Y*yY(z_);
            Zp = Zfl(z_) + X*xZ(z_) + Y*yZ(z_);
            s=surface(Xp,Yp,Zp,FIELD(:,:,iz)/max(max(max(abs(FIELD)))));
            s.EdgeColor = 'none';
            colormap(bluewhitered);
        end
    end
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
