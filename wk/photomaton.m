if 0
%% Photomaton : real space
% FIELD = ni00;   FNAME = 'ni00'; FIELDLTX = '$n_i^{00}$'; XX = RR; YY = ZZ;
% FIELD = ne00;   FNAME = 'ne00'; FIELDLTX = '$n_e^{00}$'; XX = RR; YY = ZZ;
% FIELD = dens_i; FNAME = 'ni';   FIELDLTX = '$n_i$'; XX = RR; YY = ZZ;
% FIELD = dens_e; FNAME = 'ne';   FIELDLTX = '$n_e$'; XX = RR; YY = ZZ;
% FIELD = temp_i; FNAME = 'Ti';   FIELDLTX = '$T_i$'; XX = RR; YY = ZZ;
% FIELD = temp_e; FNAME = 'Te';   FIELDLTX = '$T_e$'; XX = RR; YY = ZZ;
FIELD = phi; FNAME = 'phi'; FIELDLTX = '$\phi$'; XX = RR; YY = ZZ;
% FIELD = drphi; FNAME = 'ZF'; FIELDLTX = '$u^{ZF}_z$'; XX = RR; YY = ZZ;
% FIELD = -dr2phi; FNAME = 'shear'; FIELDLTX = '$s$'; XX = RR; YY = ZZ;
TNAME = [];
tf = 2500; [~,it1] = min(abs(Ts2D-tf)); TNAME = [TNAME,'_',num2str(tf)];
tf = 3000; [~,it2] = min(abs(Ts2D-tf)); TNAME = [TNAME,'_',num2str(tf)];
tf = 3500; [~,it3] = min(abs(Ts2D-tf)); TNAME = [TNAME,'_',num2str(tf)];
tf = 4000; [~,it4] = min(abs(Ts2D-tf)); TNAME = [TNAME,'_',num2str(tf)];
fig = figure; FIGNAME = [FNAME,TNAME,'_snaps','_',PARAMS]; set(gcf, 'Position',  [100, 100, 1500, 350])
plt = @(x) x./max(max(x));
    subplot(141)
        DATA = plt(FIELD(:,:,it1));
        pclr = pcolor((XX),(YY),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap(bluewhitered); caxis([-1,1]);
        xlabel('$r/\rho_s$'); ylabel('$z/\rho_s$');set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it1)));
    subplot(142)
        DATA = plt(FIELD(:,:,it2));
        pclr = pcolor((XX),(YY),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap(bluewhitered); caxis([-1,1]);
        xlabel('$r/\rho_s$');ylabel('$z/\rho_s$'); set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it2)));
    subplot(143)
        DATA = plt(FIELD(:,:,it3));
        pclr = pcolor((XX),(YY),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap(bluewhitered); caxis([-1,1]);
        xlabel('$r/\rho_s$');ylabel('$z/\rho_s$');set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it3)));
    subplot(144)
        DATA = plt(FIELD(:,:,it4));
        pclr = pcolor((XX),(YY),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap(bluewhitered); caxis([-1,1]);
        xlabel('$r/\rho_s$');ylabel('$z/\rho_s$'); set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it4)));
    suptitle([FIELDLTX,', $\nu_{',CONAME,'}=$', num2str(NU), ', $\eta_B=$',num2str(ETAB),...
        ', $L=',num2str(L),'$, $N=',num2str(Nr),'$, $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
        ' $\mu_{hd}=$',num2str(MU)]);
save_figure
end

if 0
%% Photomaton : quiver ExB velocity
figure
skip = 2;
plt = @(x) x./max(max(x));
FNAME = 'ZF'; FIELDLTX = '$\bm{u}^{ZF}$'; XX = RR; YY = ZZ;
tf = 1200; [~,it1] = min(abs(Ts2D-tf)); TNAME = [TNAME,'_',num2str(tf)];
UY = plt(drphi(1:skip:end,1:skip:end,it1)); UX = plt(-dzphi(1:skip:end,1:skip:end,it1)); 
pclr = pcolor(XX,YY,plt(ni00(:,:,it1))); set(pclr, 'edgecolor','none');
hold on
quiver((XX(1:skip:end,1:skip:end)),(YY(1:skip:end,1:skip:end)),UX,UY,'r'); xlim(L/2*[-1 1]); ylim(L/2*[-1 1]);
pbaspect([1 1 1])
xlabel('$r/\rho_s$');ylabel('$z/\rho_s$');set(gca,'ytick',[]); 
title(sprintf('$t c_s/R=%.0f$',Ts2D(it1)));
end


%%
if 0
%% Show frame in kspace
tf = 800; [~,it2] = min(abs(Ts2D-tf)); [~,it5] = min(abs(Ts5D-tf));
fig = figure; FIGNAME = ['krkz_',sprintf('t=%.0f',Ts2D(it2)),'_',PARAMS];set(gcf, 'Position',  [100, 100, 700, 600])
CLIM = [0,1];
    subplot(223); plt = @(x) fftshift((abs(x)),2)/max(max(abs(x)));
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(PHI(:,:,it2))); set(pclr, 'edgecolor','none');
        caxis(CLIM); colormap hot
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('$t c_s/R=%.0f$',Ts2D(it2))); legend('$|\hat\phi|$');
    subplot(222);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Ni00(:,:,it2))); set(pclr, 'edgecolor','none');
        caxis(CLIM); colormap hot
        xlabel('$k_r$'); ylabel('$k_z$'); legend('$|\hat n_i^{00}|$');
    subplot(221);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Ne00(:,:,it2))); set(pclr, 'edgecolor','none');
        caxis(CLIM); colormap hot
        xlabel('$k_r$'); ylabel('$k_z$'); legend('$|\hat n_e^{00}|$');
    subplot(224);
    colorbar;
        caxis(CLIM); colormap hot
%         pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Si00(:,:,it5))); set(pclr, 'edgecolor','none'); colorbar;
%         xlabel('$k_r$'); ylabel('$k_z$');legend('$\hat S_i^{00}$');
save_figure
end