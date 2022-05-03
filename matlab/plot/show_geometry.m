function [ FIGURE ] = show_geometry(DATA,OPTIONS)
% filtering Z pinch and torus
if DATA.Nz > 1 % Torus flux tube geometry
    Nptor  = DATA.Nz; Tgeom = 1;
    Nturns = 1;
    a      = DATA.EPS; % Torus minor radius
    R      = 1.; % Torus major radius
    q      = (DATA.Nz>1)*DATA.Q0; % Flux tube safety factor
    theta  = linspace(-pi, pi, Nptor)   ; % Poloidal angle
    phi    = linspace(0., 2.*pi, Nptor) ; % Toroidal angle
    [t, p] = meshgrid(phi, theta);
    x_tor = (R + a.*cos(p)) .* cos(t);
    y_tor = (R + a.*cos(p)) .* sin(t);
    z_tor = a.*sin(p);
    DIMENSIONS = [50 50 1200 600];
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
% xvec
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
theta  = DATA.z   ; % Poloidal angle
% Plot x,y,bvec at these points
scale = 0.10;

%% Plot plane result
OPTIONS.POLARPLOT = 0;
OPTIONS.PLAN      = 'xy';
r_o_R             = DATA.rho_o_R;
[X,Y]             = meshgrid(r_o_R*DATA.x,r_o_R*DATA.y);
max_              = 0;
FIELDS            = zeros(DATA.Ny,DATA.Nx,DATA.Nz);

FIGURE.fig = figure; FIGURE.FIGNAME = ['geometry','_',DATA.PARAMS]; set(gcf, 'Position',  DIMENSIONS)
for it_ = 1:numel(OPTIONS.TIME)
subplot(1,numel(OPTIONS.TIME),it_)
    %plot magnetic geometry
    if OPTIONS.PLT_MTOPO
    magnetic_topo=surf(x_tor, y_tor, z_tor); hold on;alpha 1.0;%light('Position',[-1 1 1],'Style','local')
    set(magnetic_topo,'edgecolor',[1 1 1]*0.7,'facecolor','none')
    end
    %plot field line
    theta  = linspace(-Nturns*pi, Nturns*pi, 512)   ; % Poloidal angle
    plot3(Xfl(theta),Yfl(theta),Zfl(theta)); hold on;
    %plot vector basis
    theta  = DATA.z   ; % Poloidal angle
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
        z_             = DATA.z(iz);    
        Xp = Xfl(z_) + X*xX(z_) + Y*yX(z_);
        Yp = Yfl(z_) + X*xY(z_) + Y*yY(z_);
        Zp = Zfl(z_) + X*xZ(z_) + Y*yZ(z_);
        s=surface(Xp,Yp,Zp,FIELDS(:,:,iz)/max_);
        s.EdgeColor = 'none';
        colormap(bluewhitered);
%         caxis([-1,1]);
    end
    if DATA.Nz == 1
        xlim([0.65 0.75]);
        ylim([-0.1 0.1]);
        zlim([-0.2 0.2]);
    end
end
    %%
    axis equal
    view([1,-2,1])

end

