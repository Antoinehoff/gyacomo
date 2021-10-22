[~,itstart] = min(abs(Ts3D-tstart));
[~,itend]   = min(abs(Ts3D-tend));
trange = itstart:itend;
%full kperp points
phi_k_2 = reshape(mean(mean(abs(FIELD(:,:,:,trange)),3),4).^2,[numel(KX),1]);
kperp  = reshape(sqrt(KX.^2+KY.^2),[numel(KX),1]);
% interpolated kperps
nk_noAA = floor(2/3*numel(kx));
kp_ip = kx;
[thg, rg] = meshgrid(linspace(0,pi,2*nk_noAA),kp_ip);
[xn,yn] = pol2cart(thg,rg);
[ky_s, sortIdx] = sort(ky);
[xc,yc] = meshgrid(ky_s,kx);
Z_rth = interp2(xc,yc,squeeze(mean((abs(FIELD(:,sortIdx,trange))).^2,3)),xn,yn);
field_kp = mean(Z_rth,2);
%for theorical trends
a1 = field_kp(2)*kp_ip(2).^(13/3);
a2 = field_kp(2)*kp_ip(2).^(3)./(1+kp_ip(2).^2).^(-2);
fig = figure; FIGNAME = ['cascade','_',FNAME,'_',PARAMS];set(gcf, 'Position',  [100, 100, 800, 300])
% scatter(kperp,phi_k_2,'.k','MarkerEdgeAlpha',0.4,'DisplayName','$|\phi_k|^2$'); hold on; grid on;
if NORMALIZED
   plt = @(x) x./max(x);
else
   plt = @(x) x;
end
plot(kp_ip,plt(field_kp),'^','DisplayName',['$\langle|',FIELDLTX,'|^2\rangle_{k_\perp}$']); hold on;
if TRENDS
plot(kp_ip,a1*kp_ip.^(-13/3),'-','DisplayName','$k^{-13/3}$');
plot(kp_ip,a2/100*kp_ip.^(-3)./(1+kp_ip.^2).^2,'-','DisplayName','$k^{-3}/(1+k^2)^2$');
end
if LOGSCALE
    set(gca, 'XScale', 'log');set(gca, 'YScale', 'log');
    xlim([0.1,10]);
end
grid on
xlabel('$k_\perp \rho_s$'); legend('show','Location','northeast')
title({['$\nu_{',CONAME,'}=$', num2str(NU), ', $\kappa_N=$',num2str(K_N),...
', $L=',num2str(L),'$, $N=',num2str(Nx),'$'];[' $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
' $\mu_{hd}=$',num2str(MU),', $',num2str(round(tstart)),'<t<',num2str(round(tend)),'$']});
save_figure
clear Z_rth phi_kp ni_kp Ti_kp
