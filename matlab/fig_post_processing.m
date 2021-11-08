%% Load figure
figpath = 'C:\Users\antoi\Desktop\gamma_eta_05_nu_1e-01_trunc';
fig = openfig(figpath);

%% Load data
axObjs = fig.Children;
dataObjs = findobj(fig,'-property','YData');
Nlines = numel(dataObjs);

%% Post processing
sigma = zeros(Nlines,1);
mu    = sigma;
tmin  = 350;
for i = 1:Nlines
    x = dataObjs(i).XData;
    
    [~, itmin] = min(abs(tmin-x));
    
    y = dataObjs(i).YData;
    
    sigma(i) = std(y(itmin:end));
    mu(i)    = mean(y(itmin:end));
end
mu = flip(mu')
sigma = flip(sigma')

%% Plot mean with error bar




if 0
%% Handwritten results for nu = 1.0, 150x75, L=100, DGGK
Results_150x75.Gamma = [0.3794    0.3194    0.3226, 0.0098    0.0221    0.0321    0.0400 0.2897    0.2886    0.2569 0.0104    0.0086    0.0276    0.0320 0.1375    0.1633    0.0848];
Results_150x75.error = [0.1909    0.1207    0.1336, 0.0028    0.0038    0.0058    0.0086 0.0832    0.0624    0.0557 0.0021    0.0023    0.0068    0.0088 0.0821    0.0278    0.0083];
Results_150x75.P     = [2,          4,          6,       2,       4,         6         8 2,          4,          6  2,          4,          6         8  2,          4,          10];
Results_150x75.J     = [1,          2,          3        1        2          3         4 1,          2,          3  1,          2,          3         4  1,          2            5];
Results_150x75.etaB  = [0.49,        0.49       0.49    0.59      0.59      0.59    0.59 0.50,     0.50       0.50  0.60,     0.60       0.60       0.60 0.51        0.51      0.51];
Results_150x75.nu    = [1.0,        1.0       1.0       1.0     1.0         1.0      1.0 0.5,        0.5       0.5  0.5,        0.5       0.5       0.5  0.1        0.1         0.1];
Results_150x75.mrkx  = [ '*',        '*',    '*',        '*',     '*',      '*',    '*', 'o',        'o',      'o', 'o',        'o',      'o',       'o' 's'        's'         's'];
Results_150x75.iclr  = [  1,        2,          3,      1          2        3          4  1,        2,          3    1,        2,          3           4  1         2             5];

% Ricci_Rogers.Gamma = [2 1e-1];
% Ricci_Rogers.etaB  = [0.5 1.0];
Ricci_Rogers.Gamma = [10  1e-2];
Ricci_Rogers.etaB  = [0.5  1.25];

if 1
% Fig 3 of Ricci Rogers 2006
SCALING = 2*sqrt(2);
fig = figure;
semilogy(Ricci_Rogers.etaB,Ricci_Rogers.Gamma,'--','color',[0,0,0]+0.6);
hold on;
plot(10,10,'color',line_colors(1,:)); plot(10,10,'color',line_colors(2,:));
plot(10,10,'color',line_colors(3,:)); plot(10,10,'color',line_colors(4,:));
plot(10,10,'color',line_colors(5,:));
plot(10,10,'*k','MarkerSize',10, 'LineWidth',1.0);
plot(10,10,'ok','MarkerSize',10, 'LineWidth',1.0);
plot(10,10,'sk','MarkerSize',10, 'LineWidth',1.0);

res = Results_150x75;
for i = 1:numel(res.Gamma)
    errorbar(res.etaB(i),res.Gamma(i)*SCALING,res.error(i)*SCALING,...
        res.mrkx(i),'DisplayName','256x128', 'color', line_colors(res.iclr(i),:),...
        'MarkerSize',12, 'LineWidth',2.0);
    hold on;
end

   xlabel('$L_n/L_B$'); ylabel('$2\sqrt(2)\Gamma^\infty_{part}$') 
end
grid on; title('$L=100$, $150\times75$, $\nu_{hyp}=0.1$'); 
xlim([0,1.75]); ylim([1e-6,100])
legend('Mix. Length, Ricci 2006','$P=2$, $J=1$','$P=4$, $J=2$','$P=6$, $J=3$','$P=8$, $J=4$','$P=10$, $J=5$',...
    '$\nu_{DGGK}=1.0$', '$\nu_{DGGK}=0.5$','$\nu_{DGGK}=0.1$');
plot([0.3 0.3],[1e-6,1e2],'r')
plot([1.6 1.6],[1e-6,1e2],'r')
plot([0.5],[0.3965],'--','color',[0,0,0]+0.6)
end