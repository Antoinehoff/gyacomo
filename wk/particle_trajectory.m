R = 1;
a = 0.05;
d = 0.02;
Wc= 80;


figure

% Magnetic field line
Bx = @(R,th) R*cos(th);
By = @(R,th) R*sin(th);
Bz = @(R,th) 0*(th);

thB = linspace(0,2*pi,128);
plot3(Bx(R,thB),By(R,thB),Bz(R,thB)); hold on;

% Particle trochoid
Tx = @(R,th) a*cos(th).*cos(Wc*th);
Ty = @(R,th) a*sin(th).*cos(Wc*th);
Tz = @(R,th) a*sin(Wc*th)+d*th;

thP = linspace(0,2.5*pi,2048);
Xp = Bx(R,thP)+Tx(R,thP);
Yp = By(R,thP)+Ty(R,thP);
Zp = Bz(R,thP)+Tz(R,thP);
plot3(Xp,Yp,Zp)
plot3(Xp(1),Yp(1),Zp(1),'xk','MarkerSize',10)
plot3(Xp(end),Yp(end),Zp(end),'ok','MarkerSize',8,'MarkerFaceColor','k')


% finish the plot
axis equal
xlim(1.2*R*[-1 1])
ylim(1.2*R*[-1 1])
zlim(0.25*[-1 1])
axis off