if 0
%% Photomaton : real space

% Chose the field to plot
% FIELD = ni00;   FNAME = 'ni00'; FIELDLTX = 'n_i^{00}';
% FIELD = ne00;   FNAME = 'ne00'; FIELDLTX = 'n_e^{00}'
% FIELD = dens_i; FNAME = 'ni';   FIELDLTX = 'n_i';
% FIELD = dens_e; FNAME = 'ne';   FIELDLTX = 'n_e';
% FIELD = temp_i; FNAME = 'Ti';   FIELDLTX = 'T_i';
% FIELD = temp_e; FNAME = 'Te';   FIELDLTX = 'T_e';
FIELD = phi; FNAME = 'phi'; FIELDLTX = '\phi';

% Chose when to plot it
tf = [0 10 50 100];

% Slice
ix = 1; iy = 1; iz = 1;
% plt = @(x,it) real(x(ix, :, :,it)); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; FIELDLTX = [FIELDLTX,'(x=',num2str(round(x(ix))),')']
% plt = @(x,it) real(x( :,iy, :,it)); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; FIELDLTX = [FIELDLTX,'(y=',num2str(round(y(iy))),')']
% plt = @(x,it) real(x( :, :,iz,it)); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; FIELDLTX = [FIELDLTX,'(x=',num2str(round(z(iz))),')'] 

% Averaged
% plt = @(x,it) mean(x(:,:,:,it),1); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; FIELDLTX = ['\langle',FIELDLTX,'\rangle_x']
% plt = @(x,it) mean(x(:,:,:,it),2); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; FIELDLTX = ['\langle',FIELDLTX,'\rangle_y']
plt = @(x,it) mean(x(:,:,:,it),3); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; FIELDLTX = ['\langle',FIELDLTX,'\rangle_z'] 


%
TNAME = [];
fig = figure; FIGNAME = [FNAME,TNAME,'_snaps','_',PARAMS]; set(gcf, 'Position',  [100, 100, 1500, 350])
plt_2 = @(x) x./max(max(x));
    for i_ = 1:numel(tf)
    [~,it] = min(abs(Ts3D-tf(i_))); TNAME = [TNAME,'_',num2str(Ts3D(it))];
    subplot(1,numel(tf),i_)
        DATA = plt_2(squeeze(plt(FIELD,it)));
        pclr = pcolor((X),(Y),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap(bluewhitered); caxis([-1,1]);
        xlabel(latexize(XNAME)); ylabel(latexize(YNAME));set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts3D(it)));
    end
    legend(latexize(FIELDLTX));
save_figure
end

if 0
%% Photomaton : quiver ExB velocity
figure
skip = 2;
plt = @(x) x./max(max(x));
FNAME = 'ZF'; FIELDLTX = '$\bm{u}^{ZF}$'; X_XY = RR; Y_XY = ZZ;
tf = 1200; [~,it1] = min(abs(Ts2D-tf)); TNAME = [TNAME,'_',num2str(tf)];
UY = plt(drphi(1:skip:end,1:skip:end,it1)); UX = plt(-dzphi(1:skip:end,1:skip:end,it1)); 
pclr = pcolor(X_XY,Y_XY,plt(ni00(:,:,it1))); set(pclr, 'edgecolor','none');
hold on
quiver((X_XY(1:skip:end,1:skip:end)),(Y_XY(1:skip:end,1:skip:end)),UX,UY,'r'); xlim(L/2*[-1 1]); ylim(L/2*[-1 1]);
pbaspect([1 1 1])
xlabel('$r/\rho_s$');ylabel('$z/\rho_s$');set(gca,'ytick',[]); 
title(sprintf('$t c_s/R=%.0f$',Ts2D(it1)));
end