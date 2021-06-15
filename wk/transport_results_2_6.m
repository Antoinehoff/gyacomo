default_plots_options
%% nuDGGK = 0.1
figure; set(gcf, 'Position',  [100, 100, 1200, 350])
% 6,3 nuDGGK=0.1
eta   = [0.5,0.6,0.7];
gamma = [116 24.396 3.684];
std_g = [23  3  0.4];
subplot(121)
errorbar(eta,gamma,std_g,'-','DisplayName','$P,J=06,03$, $\nu_{DGGK}=0.1$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(6,:)); hold on;
subplot(122)

% 10,5 nuDGGK = 0.1
eta = [0.5 0.6,0.7];
gamma = [110 24.9 4.1];
std_g = [23 4 0.38];
shear = [4.6 4.0 1.35];
std_s = [0.5 0.41 0.12];
subplot(121)
errorbar(eta,gamma,std_g,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.1$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(6,:)); hold on;
subplot(122)
errorbar(eta,shear,std_s,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.1$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(6,:)); hold on;

kmax = [0.576 0.5236 0.314];
gmax = [0.305 0.20 0.06];
subplot(121)
% Mix len Ricci 2006
semilogy([0.5  1.0],[10  1e-1],'--','color',[0,0,0]+0.6, 'DisplayName','Mix. Length');
% Bohm transport :
btransp = gmax./kmax.^2;
semilogy(.5:.1:.7,btransp','--','color',[0.4,0,0]+0.6, 'DisplayName','$\gamma/k^2$');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$\Gamma^\infty_{part}$') 
subplot(122)
semilogy([0.5 0.6 0.7],gmax,'--','color',[0,0,0]+0.6, 'DisplayName','Lin. Max Growth');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$V_E$''') 

%% nuDGGK = 0.01
figure; set(gcf, 'Position',  [100, 100, 1200, 350])
% 6,3 nuDGGK=0.01
eta = [0.5,0.6,0.7,0.8];
gamma = [30.2 4.3 0.37 0.06];
std_g = [4 0.7 0.05 0.008];
shear = [5.0 2.0 0.6 0.16];
std_s = [0.37 0.17 0.03 0.012];
subplot(121)
errorbar(eta,gamma,std_g,'-','DisplayName','$P,J=06,03$, $\nu_{DGGK}=0.01$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', 'b'); hold on;
subplot(122)
errorbar(eta,shear,std_s,'-','DisplayName','$P,J=06,03$, $\nu_{DGGK}=0.01$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', 'b'); hold on;

% % 10,5 nuDGGK=0.01
% eta    = [0.6,0.7,0.8 0.9];
% gamma  = [2.4949 0.23368 0.057792 0.023572];
% std_g  = [0.8696 0.085267 0.0060463 0.0046137];
% shear = [1.7094 0.55278 0.2054 0.085678];
% std_s = [0.2428 0.068223 0.012734 0.012291];
% subplot(121)
% errorbar(eta,gamma,std_g,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.5$',...
%     'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(5,:)); hold on;
% subplot(122)
% errorbar(eta,shear,std_s,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.5$',...
%     'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(5,:)); hold on;

kmax = [0.733 0.63 0.52 0.47];
gmax = [0.32 0.22 0.11 0.027];
subplot(121)
% Mix len Ricci 2006
semilogy([0.5  1.0],[10  1e-1],'--','color',[0,0,0]+0.6, 'DisplayName','Mix. Length');
% Bohm transport :
btransp = gmax./kmax.^2;
semilogy(.5:.1:.8,btransp','--','color',[0.4,0,0]+0.6, 'DisplayName','$\gamma/k^2$');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$\Gamma^\infty_{part}$') 
subplot(122)
semilogy([0.5:0.1:0.8],gmax,'--','color',[0,0,0]+0.6, 'DisplayName','Lin. Max Growth');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$V_E$''') 

%% nuDGGK = 0.001
figure; set(gcf, 'Position',  [100, 100, 1200, 350])
% 6,3
eta =   [0.6        0.7     0.8 0.9];
gamma = [0.26 0.088 0.042 0.0156];
std_g  = [0.4 0.06 0.003 0.0014];
gbmax  = [1.2 0.33 0.042 0.0156];
gbmin  = [0.015 0.035 0.003 0.0156];
shear = [0.75 0.55 0.34 0.19];
std_s  = [0.33 0.12 0.012];
sbmax  = [1.4 0.796 0.34 0.19];
sbmin  = [0.4 0.4 0.34 0.19];

subplot(121)
plot(eta,gamma,'-','DisplayName','$P,J=06,03$, $\nu_{DGGK}=0.001$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', 'k'); hold on;
subplot(122)
plot(eta,shear,'-','DisplayName','$P,J=06,03$, $\nu_{DGGK}=0.001$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', 'k'); hold on;

% 10,5
% eta =   [0.5        0.6        0.7      0.8       0.9];
% gamma = [12.2       0.19761 0.088 0.04253 0.026037];
% std_g   = [4.7      0.21328 0.065 0.009 0.00118];
% gbmax  = [12.2      0.8 0.25 0.04253 0.026037];
% gbmin  = [12.2      0.02 0.03 0.04253 0.026037];
% shear = [0.65361 0.46548 0.30645 0.21123 ];
% std_s  = [0.21288 0.10233 0.02886 0.0044664 ];
% sbmax  = [1 0.7 0.30645 0.21123 ];
% sbmin  = [0.4 0.35 0.30645 0.21123 ];
% errorbar(eta,gamma,std_,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.1$',...
%     'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(7,:)); hold on;
% subplot(121)
% plot(eta,gamma,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.1$',...
%     'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(7,:)); hold on;
% plot(eta,gbmax,'-.','DisplayName',' ','Color',line_colors(7,:)); hold on;
% plot(eta,gbmin,'-.','DisplayName',' ','Color',line_colors(7,:)); hold on;
% subplot(122)
% plot(eta(2:end),shear,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.1$',...
%     'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(7,:)); hold on;
% plot(eta(2:end),sbmax,'-.','DisplayName',' ','Color',line_colors(7,:)); hold on;
% plot(eta(2:end),sbmin,'-.','DisplayName',' ','Color',line_colors(7,:)); hold on;

kmax = [0.27 0.68 0.52 0.31];
gmax = [0.838 0.195 0.109 0.029];
subplot(121)
% Mix len Ricci 2006
semilogy([0.5  1.0],[10  1e-1],'--','color',[0,0,0]+0.6, 'DisplayName','Mix. Length');
% Bohm transport :
btransp = gmax./kmax.^2;
semilogy(.6:.1:.9,btransp','--','color',[0.4,0,0]+0.6, 'DisplayName','$\gamma/k^2$');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$\Gamma^\infty_{part}$') 
subplot(122)
semilogy([.6:.1:.9],gmax,'--','color',[0,0,0]+0.6, 'DisplayName','Lin. Max Growth');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$V_E$''') 
subplot(121)
plot(eta,gbmax,'-.','DisplayName',' ','Color','k'); hold on;
plot(eta,gbmin,'-.','DisplayName',' ','Color','k'); hold on;
subplot(122)
plot(eta,sbmax,'-.','DisplayName',' ','Color','k'); hold on;
plot(eta,sbmin,'-.','DisplayName',' ','Color','k'); hold on;
