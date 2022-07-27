KT_a = [3:0.5:7];
g_max= KT_a*0;
g_avg= KT_a*0;
g_std= KT_a*0;
k_max= KT_a*0;

NU    = 0.05;
DT    = 2e-3;
TMAX  = 50;
ky_   = 0.3;
SIMID = 'linear_CBC_KT_scan_ky=0.3_CLOS_0_DGGK_0.05';  % Name of the simulation
% SIMID = 'linear_CBC_PJ_scan_KT_4_ky=0.3_CLOS_0_TMAX_100';  % Name of the simulation
RUN   = 1;
figure

for P = [20]

i=1;
for K_T = KT_a
    
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
errorbar(KT_a,y_,e_,...
    'LineWidth',1.2,...
    'DisplayName',['(',num2str(P),',',num2str(P/2),')']); 
hold on;
title(['Linear CBC $K_T$ threshold $k_y=$',num2str(ky_),' (CLOS = 1)']);
legend('show'); xlabel('$K_T$'); ylabel('$\gamma$');
drawnow
end

