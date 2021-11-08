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
else % Zpinch geometry
    Nturns = 0.1; Tgeom = 0;
    a      = 1; % Torus minor radius
    R      = 0; % Torus major radius
    q      = 0; 
    theta  = linspace(-Nturns*pi, Nturns*pi, 100); % Poloidal angle
    phi    = linspace(-0.1, 0.1, 3); % cylinder height
    [t, p] = meshgrid(phi, theta);
    x_tor = a.*cos(p);
    z_tor = a.*sin(p);
    y_tor = t;
end
FIGURE.fig = figure; set(gcf, 'Position',  [50 50 1200 600])
torus=surf(x_tor, y_tor, z_tor); hold on;alpha 1.0;%light('Position',[-1 1 1],'Style','local')
set(torus,'edgecolor',[1 1 1]*0.8,'facecolor','none')
xlabel('X');ylabel('Y');zlabel('Z');
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
theta  = linspace(-Nturns*pi, Nturns*pi, 512)   ; % Poloidal angle
plot3(Xfl(theta),Yfl(theta),Zfl(theta))
% Planes plot
theta  = DATA.z   ; % Poloidal angle
plot3(Xfl(theta),Yfl(theta),Zfl(theta),'ok')
% Plot x,y,bvec at these points
scale = 0.15;
quiver3(Xfl(theta),Yfl(theta),Zfl(theta),scale*xX(theta),scale*xY(theta),scale*xZ(theta),0,'r');
quiver3(Xfl(theta),Yfl(theta),Zfl(theta),scale*yX(theta),scale*yY(theta),scale*yZ(theta),0,'g');
quiver3(Xfl(theta),Yfl(theta),Zfl(theta),scale*bX(theta),scale*bY(theta),scale*bZ(theta),0,'b');

%% Plot plane result
OPTIONS.POLARPLOT = 0;
OPTIONS.PLAN      = 'xy';
r_o_R             = DATA.rho_o_R;
[Y,X]             = meshgrid(r_o_R*DATA.y,r_o_R*DATA.x);
max_              = 0;
FIELDS            = zeros(DATA.Nx,DATA.Ny,numel(OPTIONS.PLANES));

for it_ = 1:numel(OPTIONS.TIME)
[~,it] = min(abs(OPTIONS.TIME(it_)-DATA.Ts3D));
for iz = OPTIONS.PLANES
    OPTIONS.COMP   = iz;
    toplot         = process_field(DATA,OPTIONS);
    FIELDS(:,:,iz) = toplot.FIELD(:,:,it); 
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
    caxis([-1,1]);
end
end
%%
axis equal
view([1,-2,1])

end

