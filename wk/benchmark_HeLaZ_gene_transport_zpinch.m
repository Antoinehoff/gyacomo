%% Nu = 0
KN                 = [   1.6    2.0    2.5];
Ginf_200x32x05x03  = [2.8e-3 3.1e-1 2.4e+1];
Gstd_200x32x05x03  = [4.4e-3 2.9e-1 2.0e+0];
Gmin_200x32x05x03  = [6.1e-5 2.5e-3 1.9e+1];
Gmax_200x32x05x03  = [2.3e-2 1.2e+0 2.9e+1];

Ginf_200x32x11x06  = [7.3e-3 6.3e-1 2.4e+1];
Gstd_200x32x11x06  = [1.3e-2 4.0e-1 5.9e+0];
Gmin_200x32x11x06  = [1.8e-5 4.8e-3 1.7e+1];
Gmax_200x32x11x06  = [1.7e-2 2.0e+0 3.1e+1];

Ginf_200x32x21x11  = [1.2e-2 0 0];

Ginf_GENE          = [8.5e-3 5.5e-1 2.3e+1];
Gstd_GENE          = [5.1e-3 3.7e-1 2.6e+0];
Gmin_GENE          = [8.0e-4 9.2e-3 2.0e+1];
Gmax_GENE          = [6.0e-3 1.5e+0 2.3e+1];


fig = figure; set(gcf,'Position',[250 250 600 300]);
% plot(KN,Ginf_200x32x05x03./Ginf_GENE,'v','DisplayName','HeLaZ 200x32x05x03','MarkerSize',10); hold on;
% plot(KN,Ginf_200x32x11x06./Ginf_GENE,'^','DisplayName','HeLaZ 200x32x10x05','MarkerSize',10); hold on;
% errorbar(KN-0.05,Ginf_200x32x05x03,Gstd_200x32x11x06/2,'v','DisplayName','HeLaZ 200x32x05x03','MarkerSize',10); hold on;
% errorbar(KN+0.05,Ginf_200x32x11x06,Gstd_200x32x11x06/2,'v','DisplayName','HeLaZ 200x32x05x03','MarkerSize',10); hold on;
% errorbar(KN     ,Ginf_GENE        ,Gstd_GENE/2        ,'xk','DisplayName','HeLaZ 200x32x05x03','MarkerSize',10); hold on;
% errorbar(KN-0.05,Ginf_200x32x05x03,Gmin_200x32x05x03,Gmax_200x32x05x03,'v','DisplayName','HeLaZ 200x32x05x03','MarkerSize',10); hold on;
% errorbar(KN+0.05,Ginf_200x32x11x06,Gmin_200x32x11x06,Gmax_200x32x11x06,'v','DisplayName','HeLaZ 200x32x05x03','MarkerSize',10); hold on;
% errorbar(KN     ,Ginf_GENE        ,Gmin_GENE        ,Gmax_GENE        ,'xk','DisplayName','HeLaZ 200x32x05x03','MarkerSize',10); hold on;
% xlabel('$\kappa_N$');
% ylabel('$\Gamma_x^\infty$');

loglog(Ginf_GENE,Ginf_200x32x05x03,'*','DisplayName','HeLaZ 5x3','MarkerSize',10); hold on;
loglog(Ginf_GENE,Ginf_200x32x11x06,'*','DisplayName','HeLaZ 10x5','MarkerSize',8); hold on;
loglog(Ginf_GENE,Ginf_200x32x21x11,'*','DisplayName','HeLaZ 20x10','MarkerSize',8); hold on;
% title('Collisionless transport benchmark');
xy = [1e-3 1e+2];
loglog(xy,xy,'--k','DisplayName','x=y');
% loglog(Ginf_GENE,Gmin_GENE,'-k','DisplayName','x=y');
% loglog(Ginf_GENE,Gmax_GENE,'-k','DisplayName','x=y');
xlabel('Gene $\Gamma_x^\infty$');
ylabel('HeLaZ $\Gamma_x^\infty$');
xlim(xy); ylim(xy);
grid on
% saveas(fig,'/home/ahoffman/Dropbox/Applications/Overleaf/Paper 1/figures/transport_benchmark.eps')


%% Nu = 0.1
fig = figure; set(gcf,'Position',[250 250 600 300]);

KN                 = [   1.6    2.0    2.5];
LDGK_200x32x05x03  = [2.6e-1 4.3e+0 3.8e+1];
LRGK_200x32x05x03  = [4.0e-1 2.4e+0 2.9e+1];
SGGK_200x32x05x03  = [1.4e-2 1.4e+0 3.5e+1];
DGGK_200x32x05x03  = [3.8e-3 1.8e-1 3.0e+1];

SGGK_200x32x11x06  = [1.5e-2 1.3e+0 0.0e+0];
DGGK_200x32x11x06  = [2.4e-3 3.8e-1 0.0e+0];

X = Ginf_GENE;
% X = KN;
clr_ = line_colors;

loglog(X,DGGK_200x32x05x03,'v','DisplayName','Dougherty','MarkerSize',10,'Color',clr_(2,:)); hold on;
loglog(X,SGGK_200x32x05x03,'s','DisplayName',   'Sugama','MarkerSize',10,'Color',clr_(1,:)); hold on;
loglog(X,LDGK_200x32x05x03,'d','DisplayName',   'Landau','MarkerSize',10,'Color',clr_(5,:)); hold on;
loglog(X,LRGK_200x32x05x03,'^','DisplayName',  'Lorentz','MarkerSize',10,'Color',clr_(3,:)); hold on;

loglog(X,DGGK_200x32x11x06,'vk','DisplayName',   'conv','MarkerSize',8); hold on;
loglog(X,SGGK_200x32x11x06,'sk','DisplayName',   'conv','MarkerSize',8); hold on;

%% Nu = 0.1
fig = figure; set(gcf,'Position',[250 250 600 300]);

% KN                 = [   1.6    2.0    2.5];
% LDGK_200x32x05x03  = [2.6e-1 3.6e+0 3.8e+1];
% LRGK_200x32x05x03  = [4.0e-1 2.0e+0 2.9e+1];
% SGGK_200x32x05x03  = [1.5e-2 1.4e+0 3.5e+1];
% DGGK_200x32x05x03  = [3.6e-3 1.8e-1 3.2e+1];
% NOCO_200x32x05x03  = [3.2e-3 3.1e-1 2.4e+1];

% LRGK_200x32x09x05  = [0.0e+0 2.0e+0 0.0e+0];
% SGGK_200x32x11x06  = [1.5e-2 1.3e+0 0.0e+0];
% DGGK_200x32x11x06  = [2.4e-3 3.8e-1 0.0e+0];

KN                = [   1.5    1.6    1.7    1.8    1.9    2.0    2.1    2.2    2.3    2.4    2.5];
LDGK_200x32x05x03 = [2.6e-2 2.6e-1 6.2e-1 1.2e+0 2.3e+0 3.6e+0 7.4e+0 1.1e+1 1.5e+1 2.3e+1 3.8e+1]; % validate with two box sizes (120 and 240)
LRGK_200x32x05x03 = [0.0e+0 4.0e-1 7.2e-1 1.2e+0 1.7e+0 2.3e+0 4.0e+0 6.0e+0 1.1e+1 1.9e+1 2.9e+1];
SGGK_200x32x05x03 = [0.0e-9 1.5e-2 9.5e-2 3.5e-1 8.0e-1 1.4e+0 2.5e+0 4.2e+0 8.0e+0 1.6e+1 3.5e+1];
hacked_sugama     = [0.0e+0 7.0e-1 0.0e+0 2.5e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0];
DGGK_200x32x05x03 = [2.0e-4 5.0e-3 9.2e-3 3.5e-2 1.3e-1 3.5e-1 6.0e-1 6.0e+0 1.1e+1 1.3e+1 3.2e+1]; % validate with two box sizes (120 and 240)
LBGK_200x32x05x03 = [0.0e+0 0.0e-1 0.0e+0 0.0e+0 0.0e+0 2.8e-1 0.0e+0 0.0e+0 0.0e+0 0.0e+0 0.0e+0]; % validate with two box sizes (120 and 240)
NOCO_200x32x05x03 = [0.0e+0 1.2e-2 4.0e-2 9.0e-2 2.3e-1 3.1e-1 7.0e-1 3.9e+0 6.9e+0 1.2e+1 2.4e+1];


% X = Ginf_GENE;
X = KN;
clr_ = line_colors;
plt = @(x) x./X;
% plt = @(x) x;
msize = 10;
semilogy(X,plt(DGGK_200x32x05x03),'v','DisplayName','Dougherty','MarkerSize',msize,'Color',clr_(2,:)); hold on;
semilogy(X,plt(SGGK_200x32x05x03),'s','DisplayName',   'Sugama','MarkerSize',msize,'Color',clr_(1,:)); hold on;
semilogy(X,plt(LDGK_200x32x05x03),'d','DisplayName',   'Landau','MarkerSize',msize,'Color',clr_(5,:)); hold on;
semilogy(X,plt(LRGK_200x32x05x03),'^','DisplayName',  'Lorentz','MarkerSize',msize,'Color',clr_(3,:)); hold on;
semilogy(X,plt(hacked_sugama),'s','DisplayName',  'hacked Sugama','MarkerSize',msize,'Color',clr_(6,:)); hold on;
semilogy(X,plt(NOCO_200x32x05x03),'*k','DisplayName',  '$\nu = 0$','MarkerSize',msize); hold on;

X = [   1.6    2.0    2.5];
% loglog(X,DGGK_200x32x11x06./X,'vk','DisplayName',   'conv','MarkerSize',8); hold on;
% loglog(X,SGGK_200x32x11x06./X,'sk','DisplayName',   'conv','MarkerSize',8); hold on;
% loglog(X,LRGK_200x32x09x05./X,'^k','DisplayName',   'conv','MarkerSize',8); hold on;
ylabel('$\Gamma_x^\infty/\kappa_N$, $\nu=0.1$');
xlabel('$\kappa_N$'); grid on
xlim([ 1.55 2.55]);
%%
% saveas(fig,'/home/ahoffman/Dropbox/Applications/Overleaf/Paper 1/figures/coll_transport_benchmark.eps')

%% Burst study Kn = 1.7
Kn = 1.7;
NU                = [0.0e+0 1.0e-2 2.0e-2 3.0e-2 4.0e-2 5.0e-2 6.0e-2 7.0e-2 8.0e-2 9.0e-2 1.0e-1];
SGGK_200x32x05x03 = [1.0e-2 1.1e-2 1.7e-2 2.3e-2 3.1e-2 4.1e-2 5.1e-2 6.1e-2 8.3e-2 8.8e-2 1.0e-1];
DGGK_200x32x05x03 = [1.0e-2 1.0e-3 1.1e-3 2.3e-3 3.0e-3 3.0e-3 3.7e-3 4.0e-3 5.7e-3 5.8e-3 9.2e-3];
LDGK_200x32x05x03 = [1.0e-2 9.0e-2 1.5e-1 2.2e-1 2.8e-1 3.5e-1 4.0e-1 4.5e-1 5.0e-1 5.6e-1 6.2e-1];
% NOCO_200x32x05x03 = [0.0e+0 0.0e+0 0.0e+0 0.0e+0];
figure
plot(NU,SGGK_200x32x05x03/Kn,'s','Color',clr_(1,:)); hold on
plot(NU,DGGK_200x32x05x03/Kn,'v','Color',clr_(2,:)); hold on
plot(NU,LDGK_200x32x05x03/Kn,'d','Color',clr_(5,:)); hold on
grid on; xlabel('$\nu$'); ylabel('$\Gamma_x^\infty/\kappa_N$');
