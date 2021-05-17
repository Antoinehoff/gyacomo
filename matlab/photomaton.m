if 0
%% Photomaton : real space
% FIELD = ni00; FNAME = 'ni'; XX = RR; YY = ZZ;
% FIELD = phi; FNAME = 'phi'; XX = RR; YY = ZZ;
FIELD = dr2phi; FNAME = 'dr2phi'; XX = RR; YY = ZZ;
% FIELD = fftshift(abs(Ni00),2); FNAME = 'Fni'; XX = fftshift(KR,2); YY = fftshift(KZ,2);
% FIELD = fftshift(abs(PHI),2);  FNAME = 'Fphi'; XX = fftshift(KR,2); YY = fftshift(KZ,2);
tf = 150;  [~,it1] = min(abs(Ts2D-tf));
tf = 300;  [~,it2] = min(abs(Ts2D-tf)); 
tf = 800; [~,it3] = min(abs(Ts2D-tf));
tf = 1200; [~,it4] = min(abs(Ts2D-tf));
fig = figure; FIGNAME = [FNAME,'_snaps','_',PARAMS]; set(gcf, 'Position',  [100, 100, 1500, 350])
plt = @(x) x;%./max(max(x));
    subplot(141)
        DATA = plt(FIELD(:,:,it1));
        pclr = pcolor((XX),(YY),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap gray
        xlabel('$r/\rho_s$'); ylabel('$z/\rho_s$');set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it1)));
    subplot(142)
        DATA = plt(FIELD(:,:,it2));
        pclr = pcolor((XX),(YY),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap gray
        xlabel('$r/\rho_s$');ylabel('$z/\rho_s$'); set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it2)));
    subplot(143)
        DATA = plt(FIELD(:,:,it3));
        pclr = pcolor((XX),(YY),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap gray
        xlabel('$r/\rho_s$');ylabel('$z/\rho_s$');set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it3)));
    subplot(144)
        DATA = plt(FIELD(:,:,it4));
        pclr = pcolor((XX),(YY),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap gray
        xlabel('$r/\rho_s$');ylabel('$z/\rho_s$'); set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it4)));
suptitle(['$\',FNAME,'$, $\nu_{',CONAME,'}=$', num2str(NU), ', $\eta_B=$',num2str(ETAB),...
    ', $P=',num2str(PMAXI),'$, $J=',num2str(JMAXI),'$']);
save_figure
end

%%
if 0
%% Show frame in kspace
tf = 1000; [~,it2] = min(abs(Ts2D-tf)); [~,it5] = min(abs(Ts5D-tf));
fig = figure; FIGNAME = ['krkz_',sprintf('t=%.0f',Ts2D(it2)),'_',PARAMS];set(gcf, 'Position',  [100, 100, 700, 600])
CLIM = [0,1e3]
    subplot(223); plt = @(x) fftshift((abs(x)),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(PHI(:,:,it2))); set(pclr, 'edgecolor','none');
        caxis(CLIM); colormap hot
        xlabel('$k_r$'); ylabel('$k_z$'); title(sprintf('$t c_s/R=%.0f$',Ts2D(it2))); legend('$|\hat\phi|$');
    subplot(222); plt = @(x) fftshift(abs(x),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Ni00(:,:,it2))); set(pclr, 'edgecolor','none');
        caxis(CLIM); colormap hot
        xlabel('$k_r$'); ylabel('$k_z$'); legend('$|\hat n_i^{00}|$');
    subplot(221); plt = @(x) fftshift(abs(x),2);
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Ne00(:,:,it2))); set(pclr, 'edgecolor','none');
        caxis(CLIM); colormap hot
        xlabel('$k_r$'); ylabel('$k_z$'); legend('$|\hat n_e^{00}|$');
    subplot(224); plt = @(x) fftshift((abs(x)),2);
    colorbar;
        caxis(CLIM); colormap hot
%         pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),plt(Si00(:,:,it5))); set(pclr, 'edgecolor','none'); colorbar;
%         xlabel('$k_r$'); ylabel('$k_z$');legend('$\hat S_i^{00}$');
save_figure
end