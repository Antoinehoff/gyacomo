% trange = 100:200;
figure;
plt = @(x) squeeze(mean(mean(real(x(:,:,trange)),2),3))/max(abs(squeeze(mean(mean(real(x(:,:,trange)),2),3))));
plot(r,plt(-dr2phi),'-k','DisplayName','Zonal shear'); hold on;
plot(r,plt(-drphi.*dzphi),'-','DisplayName','$\Pi_\phi$'); hold on;
% plot(r,plt(-drphi.*dzTe),'-','DisplayName','$\Pi_{Te}$'); hold on;
plot(r,plt(-drphi.*dzTi),'-','DisplayName','$\Pi_{Ti}$'); hold on;
plot(r,plt(-drphi.*dzphi-drphi.*dzTi),'-','DisplayName','$\Pi_\phi+\Pi_{Ti}$'); hold on;
% plot(r,plt(-drphi.*dzphi-drphi.*dzTi-drphi.*dzTe),'-','DisplayName','$\Pi_\phi+\Pi_{Te}+\Pi_{Ti}$'); hold on;
xlim([-L/2,L/2]); xlabel('$x/\rho_s$'); grid on; legend('show')