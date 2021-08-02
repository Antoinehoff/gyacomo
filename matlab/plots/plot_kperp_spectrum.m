[~,itstart] = min(abs(Ts3D-tstart));
[~,itend]   = min(abs(Ts3D-tend));
trange = itstart:itend;
%full kperp points
phi_k_2 = reshape(mean(mean(abs(PHI(:,:,:,trange)),3),4).^2,[numel(KX),1]);
kperp  = reshape(sqrt(KX.^2+KY.^2),[numel(KX),1]);
% interpolated kperps
nk_noAA = floor(2/3*numel(kx));
kp_ip = kx;
[thg, rg] = meshgrid(linspace(0,pi,2*nk_noAA),kp_ip);
[xn,yn] = pol2cart(thg,rg);
[ky_s, sortIdx] = sort(ky);
[xc,yc] = meshgrid(ky_s,kx);
Z_rth = interp2(xc,yc,squeeze(mean((abs(PHI(:,sortIdx,trange))).^2,3)),xn,yn);
phi_kp = mean(Z_rth,2);
Z_rth = interp2(xc,yc,squeeze(mean((abs(Ni00(:,sortIdx,trange))).^2,3)),xn,yn);
ni00_kp = mean(Z_rth,2);
Z_rth = interp2(xc,yc,squeeze(mean((abs(Ne00(:,sortIdx,trange))).^2,3)),xn,yn);
ne00_kp = mean(Z_rth,2);
%for theorical trends
a1 = phi_kp(2)*kp_ip(2).^(13/3);
a2 = phi_kp(2)*kp_ip(2).^(3)./(1+kp_ip(2).^2).^(-2);
fig = figure; FIGNAME = ['cascade','_',PARAMS];set(gcf, 'Position',  [100, 100, 500, 500])
% scatter(kperp,phi_k_2,'.k','MarkerEdgeAlpha',0.4,'DisplayName','$|\phi_k|^2$'); hold on; grid on;
plot(kp_ip,phi_kp,'^','DisplayName','$\langle|\phi_k|^2\rangle_{k_\perp}$'); hold on;
plot(kp_ip,ni00_kp, '^','DisplayName','$\langle|N_{i,k}^{00}|^2\rangle_{k_\perp}$'); hold on;
plot(kp_ip,ne00_kp, '^','DisplayName','$\langle|N_{e,k}^{00}|^2\rangle_{k_\perp}$'); hold on;
plot(kp_ip,a1*kp_ip.^(-13/3),'-','DisplayName','$k^{-13/3}$');
plot(kp_ip,a2/100*kp_ip.^(-3)./(1+kp_ip.^2).^2,'-','DisplayName','$k^{-3}/(1+k^2)^2$');
set(gca, 'XScale', 'log');set(gca, 'YScale', 'log'); grid on
xlabel('$k_\perp \rho_s$'); legend('show','Location','southwest')
title({['$\nu_{',CONAME,'}=$', num2str(NU), ', $\eta=$',num2str(ETAB/ETAN),...
', $L=',num2str(L),'$, $N=',num2str(Nx),'$'];[' $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
' $\mu_{hd}=$',num2str(MU)]});
xlim([0.1,10]);
save_figure
clear Z_rth phi_kp ni_kp Ti_kp