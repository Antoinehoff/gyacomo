function [ FIGURE ] = show_geometry(DATA,OPTIONS)
% filtering Z pinch and torus
if DATA.grids.Nz > 1 % Torus flux tube geometry
    Nturns = floor(abs(DATA.grids.z(1)/pi));
    Nptor  = ceil(DATA.grids.Nz*2/Nturns); Tgeom = 1;
    a      = DATA.inputs.EPS; % Torus minor radius
    R      = 1.; % Torus major radius
    q      = (DATA.grids.Nz>1)*DATA.inputs.Q0; % Flux tube safety factor
    theta  = linspace(-pi, pi, Nptor)   ; % Poloidal angle
    phi    = linspace(0., 2.*pi, Nptor) ; % Toroidal angle
    [t, p] = meshgrid(phi, theta);
    x_tor = (R + a.*cos(p)) .* cos(t);
    y_tor = (R + a.*cos(p)) .* sin(t);
    z_tor = a.*sin(p);
    DIMENSIONS = [600 600 1200 600];
else % Zpinch geometry
    Nturns = 0.1; Tgeom = 0;
    a      = 0.7; % Torus minor radius
    R      = 0; % Torus major radius
    q      = 0; 
    theta  = linspace(-Nturns*pi, Nturns*pi, 100); % Poloidal angle
    phi    = linspace(-0.1, 0.1, 5); % cylinder height
    [t, p] = meshgrid(phi, theta);
    x_tor = a.*cos(p);
    z_tor = a.*sin(p);
    y_tor = t;
    DIMENSIONS = [50 50 numel(OPTIONS.TIME)*400 500];
end
% field line coordinates
Xfl = @(z) (R+a*cos(z)).*cos(q*z);
Yfl = @(z) (R+a*cos(z)).*sin(q*z);
Zfl = @(z) a*sin(z);
Rvec= @(z) [Xfl(z); Yfl(z); Zfl(z)];
% xvec shearless
xX  = @(z) (Xfl(z)-R*cos(q*z))./sqrt((Xfl(z)-R*cos(q*z)).^2+(Yfl(z)-R*sin(q*z)).^2+Zfl(z).^2);
xY  = @(z) (Yfl(z)-R*sin(q*z))./sqrt((Xfl(z)-R*cos(q*z)).^2+(Yfl(z)-R*sin(q*z)).^2+Zfl(z).^2);
xZ  = @(z)              Zfl(z)./sqrt((Xfl(z)-R*cos(q*z)).^2+(Yfl(z)-R*sin(q*z)).^2+Zfl(z).^2);
xvec= @(z) [xX(z); xY(z); xZ(z)];
% bvec
bX  = @(z) Tgeom*(a*cos(z).*cos(q*z) - q*Yfl(z))./sqrt(Tgeom*(a*cos(z).*cos(q*z) - q*Yfl(z)).^2+(a*cos(z).*sin(q*z) + q*Xfl(z)).^2+(a*cos(z)).^2);
bY  = @(z)       (a*cos(z).*sin(q*z) + q*Xfl(z))./sqrt(Tgeom*(a*cos(z).*cos(q*z) - q*Yfl(z)).^2+(a*cos(z).*sin(q*z) + q*Xfl(z)).^2+(a*cos(z)).^2);
bZ  = @(z)                              a*cos(z)./sqrt(Tgeom*(a*cos(z).*cos(q*z) - q*Yfl(z)).^2+(a*cos(z).*sin(q*z) + q*Xfl(z)).^2+(a*cos(z)).^2);
bvec= @(z) [bX(z); bY(z); bZ(z)];
% yvec = b times x
yX  = @(z) bY(z).*xZ(z)-bZ(z).*xY(z)./sqrt((bY(z).*xZ(z)-bZ(z).*xY(z)).^2+(bZ(z).*xX(z)-bX(z).*xZ(z)).^2+(bX(z).*xY(z)-bY(z).*xX(z)).^2);
yY  = @(z) bZ(z).*xX(z)-bX(z).*xZ(z)./sqrt((bY(z).*xZ(z)-bZ(z).*xY(z)).^2+(bZ(z).*xX(z)-bX(z).*xZ(z)).^2+(bX(z).*xY(z)-bY(z).*xX(z)).^2);
yZ  = @(z) bX(z).*xY(z)-bY(z).*xX(z)./sqrt((bY(z).*xZ(z)-bZ(z).*xY(z)).^2+(bZ(z).*xX(z)-bX(z).*xZ(z)).^2+(bX(z).*xY(z)-bY(z).*xX(z)).^2);
yvec= @(z) [yX(z); yY(z); yZ(z)];
% Plot high res field line
% Planes plot
theta  = DATA.grids.z   ; % Poloidal angle
% Plot x,y,bvec at these points
scale = 0.10;

%% Plot plane result
OPTIONS.POLARPLOT = 0;
OPTIONS.PLAN      = 'xy';
r_o_R             = DATA.rho_o_R;
[X,Y]             = meshgrid(r_o_R*DATA.grids.x,r_o_R*DATA.grids.y);
max_              = 0;
FIELDS            = zeros(DATA.grids.Ny,DATA.grids.Nx,DATA.grids.Nz);

FIGURE.fig = figure; FIGURE.FIGNAME = ['geometry','_',DATA.params_string]; set(gcf, 'Position',  DIMENSIONS)
for it_ = 1:numel(OPTIONS.TIME)
subplot(1,numel(OPTIONS.TIME),it_)
    %plot magnetic geometry
    if OPTIONS.PLT_MTOPO
    magnetic_topo=surf(x_tor, y_tor, z_tor); hold on;alpha 0.5;%light('Position',[-1 1 1],'Style','local')
    set(magnetic_topo,'edgecolor',[1     1 1]*0.8,'facecolor','none')
%     set(magnetic_topo,'edgecolor','none','facecolor','white')
    end
    %plot field line
    theta  = linspace(-Nturns*pi, Nturns*pi, 512)   ; % Poloidal angle
    plot3(Xfl(theta),Yfl(theta),Zfl(theta)); hold on;
    %plot fluxe tube
    if OPTIONS.PLT_FTUBE
    theta  = linspace(-Nturns*pi, Nturns*pi, 64)    ; % Poloidal angle
    %store the shifts in an order (top left to bottom right)
    s_x = r_o_R*[DATA.x(1) DATA.x(end) DATA.x(1)   DATA.x(end)]; 
    s_y = r_o_R*[DATA.y(1) DATA.y(1)   DATA.y(end) DATA.y(end)];
    for i_ = 1:4
    vx_ = Xfl(theta) + s_x(i_)*xX(theta) + s_y(i_)*yX(theta);
    vy_ = Yfl(theta) + s_x(i_)*xY(theta) + s_y(i_)*yY(theta);
    vz_ = Zfl(theta) + s_x(i_)*xZ(theta) + s_y(i_)*yZ(theta);
    plot3(vx_,vy_,vz_,'-','color',[1.0 0.6 0.6]*0.8,'linewidth',1.5); hold on;
    end
    end
    %plot vector basis
    theta  = DATA.grids.z   ; % Poloidal angle
    plot3(Xfl(theta),Yfl(theta),Zfl(theta),'ok'); hold on;
    quiver3(Xfl(theta),Yfl(theta),Zfl(theta),scale*xX(theta),scale*xY(theta),scale*xZ(theta),0,'r');
    quiver3(Xfl(theta),Yfl(theta),Zfl(theta),scale*yX(theta),scale*yY(theta),scale*yZ(theta),0,'g');
    quiver3(Xfl(theta),Yfl(theta),Zfl(theta),scale*bX(theta),scale*bY(theta),scale*bZ(theta),0,'b');
    xlabel('X');ylabel('Y');zlabel('Z');
    %Plot time dependent data
    for iz = OPTIONS.PLANES
        OPTIONS.COMP   = iz;
        toplot         = process_field(DATA,OPTIONS);
        FIELDS(:,:,iz) = toplot.FIELD(:,:,it_); 
        tmp            = max(max(abs(FIELDS(:,:,iz))));
        if (tmp > max_) max_ = tmp;
    end
    for iz = OPTIONS.PLANES
        z_             = DATA.grids.z(iz);    
        Xp = Xfl(z_) + X*xX(z_) + Y*yX(z_);
        Yp = Yfl(z_) + X*xY(z_) + Y*yY(z_);
        Zp = Zfl(z_) + X*xZ(z_) + Y*yZ(z_);
        s=surface(Xp,Yp,Zp,FIELDS(:,:,iz)/max_);
        s.EdgeColor = 'none';
        colormap(bluewhitered);
%         caxis([-1,1]);
    end
    if DATA.grids.Nz == 1
        xlim([0.65 0.75]);
        ylim([-0.1 0.1]);
        zlim([-0.2 0.2]);
    end
end
    %%
    axis equal
    view([1,-2,1])
    grid on
end

