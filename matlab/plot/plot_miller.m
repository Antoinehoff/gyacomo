function [R,Z] = plot_miller(geom,rho,Nz,PLOT)
%
R0 = geom.R0; Z0 = geom.Z0; 
kappa = geom.kappa; delta = geom.delta; zeta = geom.zeta;

theta = linspace(0,2*pi,Nz+1);
R = R0 + rho*cos(theta+asin(delta*sin(theta)));
Z = Z0 + kappa*rho*sin(theta + zeta*sin(2*theta));

if PLOT
figure
plot(R,Z,'-b');
xlabel('R [m]');
ylabel('Z [m]');
axis tight
axis equal
end
end