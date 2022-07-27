% NU_a = [0.05 0.15 0.25 0.35 0.45];
NU_a = [0:0.05:0.5];
g_max= NU_a*0;
g_avg= NU_a*0;
g_std= NU_a*0;
k_max= NU_a*0;
CO      = 'LR';

K_T   = 5.3;
DT    = 2e-3;
TMAX  = 30;
ky_   = 0.3;
SIMID = 'linear_CBC_nu_scan_ky=0.3_CLOS_0_LRGK';  % Name of the simulation
RUN   = 1;
figure

for P = [2 4 6 10 12 20]

i=1;
for NU = NU_a
    
    quick_run
    
    g_max(i) = gmax;
    k_max(i) = kmax;
    
    g_avg(i) = lg.avg_g;
    g_std(i) = lg.std_g;
    
    i = i + 1;
end
%%

% plot(KT_a,max(g_max,0));
y_ = g_avg; 
e_ = g_std;

y_ = y_.*(y_-e_>0);
e_ = e_ .* (y_>0);
errorbar(NU_a,y_,e_,...
    'LineWidth',1.2,...
    'DisplayName',['(',num2str(P),',',num2str(P/2),')']); 
hold on;
title(['$\kappa_T=$',num2str(K_T),' $k_y=$',num2str(ky_),' (CLOS = 0)']);
legend('show'); xlabel('$\nu_{DGGK}$'); ylabel('$\gamma$');
drawnow
end

