if 0
%% Photomaton : real space
% FIELD = ni00; FNAME = 'ni'; FIELDLTX = '$N_i^{00}$'; XX = RR; YY = ZZ;
% FIELD = ne00; FNAME = 'ne'; FIELDLTX = '$N_e^{00}$'; XX = RR; YY = ZZ;
FIELD = phi; FNAME = 'phi'; FIELDLTX = '$\phi$'; XX = RR; YY = ZZ;
% FIELD = dr2phi; FNAME = 'dr2phi'; XX = RR; YY = ZZ;
tf = 100;  [~,it1] = min(abs(Ts2D-tf));
tf = 285;  [~,it2] = min(abs(Ts2D-tf)); 
tf = 350; [~,it3] = min(abs(Ts2D-tf));
tf = 1100; [~,it4] = min(abs(Ts2D-tf));
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
    suptitle([FIELDLTX,', $\nu_{',CONAME,'}=$', num2str(NU), ', $\eta_B=$',num2str(ETAB),...
        ', $L=',num2str(L),'$, $N=',num2str(Nr),'$, $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
        ' $\mu_{hd}=$',num2str(MU)]);
save_figure
end

if 0
%% Photomaton : real space
% FIELD = Ni00; FNAME = 'fni00'; FIELDLTX = '$\tilde N_i^{00}$'; XX = RR; YY = ZZ;
% FIELD = Ne00; FNAME = 'fne00'; FIELDLTX = '$\tilde N_e^{00}$'; XX = RR; YY = ZZ;
FIELD = PHI; FNAME = 'fphi'; FIELDLTX = '$\tilde \phi$'; XX = RR; YY = ZZ;
tf = 700;  [~,it1] = min(abs(Ts2D-tf));
tf = 800;  [~,it2] = min(abs(Ts2D-tf)); 
tf = 900; [~,it3] = min(abs(Ts2D-tf));
tf = 1500; [~,it4] = min(abs(Ts2D-tf));
XLIM = [0,1.0]; YLIM = [-.5;.5];
fig = figure; FIGNAME = [FNAME,'_snaps','_',PARAMS]; set(gcf, 'Position',  [100, 100, 1500, 400])
plt = @(x) fftshift((abs(x)),2)/max(max(abs(x)));
    subplot(141)
        DATA = plt(FIELD(:,:,it1));
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap gray
        xlabel('$k_r \rho_s$'); ylabel('$k_z \rho_s$'); xlim(XLIM); ylim(YLIM); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it1)));
    subplot(142)
        DATA = plt(FIELD(:,:,it2));
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap gray
        xlabel('$k_r \rho_s$'); ylabel('$k_z \rho_s$'); set(gca,'ytick',[]);  xlim(XLIM); ylim(YLIM); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it2)));
    subplot(143)
        DATA = plt(FIELD(:,:,it3));
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap gray
        xlabel('$k_r \rho_s$'); ylabel('$k_z \rho_s$');set(gca,'ytick',[]);  xlim(XLIM); ylim(YLIM); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it3)));
    subplot(144)
        DATA = plt(FIELD(:,:,it4));
        pclr = pcolor(fftshift(KR,2),fftshift(KZ,2),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap gray
        xlabel('$k_r \rho_s$'); ylabel('$k_z \rho_s$'); set(gca,'ytick',[]);  xlim(XLIM); ylim(YLIM); 
        title(sprintf('$t c_s/R=%.0f$',Ts2D(it4)));
    suptitle([FIELDLTX,', $\nu_{',CONAME,'}=$', num2str(NU), ', $\eta_B=$',num2str(ETAB),...
        ', $L=',num2str(L),'$, $N=',num2str(Nr),'$, $(P,J)=(',num2str(PMAXI),',',num2str(JMAXI),')$,',...
        ' $\mu_{hd}=$',num2str(MU)]);
save_figure
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