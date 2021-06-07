default_plots_options
%% nuDGGK = 0.1
figure; set(gcf, 'Position',  [100, 100, 1200, 350])
% 6,3 nuDGGK=1.0
eta   = [0.5,0.6,0.7,0.8];
gamma = [32.6443 3.6895 0.3744 0.056];
std_g = [6.1529 0.7986 0.0493 0.0134];
subplot(121)
errorbar(eta,gamma,std_g,'-','DisplayName','$P,J=06,03$, $\nu_{DGGK}=1.0$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(1,:)); hold on;
subplot(122)

% 10,5 nuDGGK = 1.0
eta = [0.5 0.6,0.7,0.8];
gamma = [32.6 3.5546 0.3917 0.0500];
std_g = [7.7 0.5846 0.0486 0.0088];
shear = [1.8505 0.60866 0.048249];
std_s = [0.1599 0.00614 0.0061403];
subplot(121)
errorbar(eta,gamma,std_g,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=1.0$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(2,:)); hold on;
subplot(122)
errorbar(eta(2:end),shear,std_s,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=1.0$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(2,:)); hold on;

% 12,6 nuDGGK = 1.0
eta = [0.6];
gamma = [4.064];
std_g = [0.7964];
subplot(121)
errorbar(eta,gamma,std_g,'-','DisplayName','$P,J=12,06$, $\nu_{DGGK}=1.0$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(3,:)); hold on;
subplot(122)


% Mix len Ricci 2006
subplot(121)
semilogy([0.5  1.0],[10  1e-1],'--','color',[0,0,0]+0.6, 'DisplayName','Mix. Length');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$\Gamma^\infty_{part}$') 
subplot(122)
grate = [0.2131 0.106 0.02021];
semilogy([0.6 0.7 0.8],grate,'--','color',[0,0,0]+0.6, 'DisplayName','Lin. Max Growth');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$V_E$''') 

%% nuDGGK = 0.5
figure; set(gcf, 'Position',  [100, 100, 1200, 350])
% 6,3 nuDGGK=0.5
eta = [0.5,0.6,0.7,0.8];
gamma = [20.511 2.6292 0.2353 0.057];
std_g = [3.67 1.2 0.055 0.004];
shear = [nan 1.7417  0.57345 0.25155];
std_s = [nan 0.35991 0.041 0.00913];
subplot(121)
errorbar(eta,gamma,std_g,'-','DisplayName','$P,J=06,03$, $\nu_{DGGK}=0.5$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(4,:)); hold on;
subplot(122)
errorbar(eta,shear,std_s,'-','DisplayName','$P,J=06,03$, $\nu_{DGGK}=0.5$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(4,:)); hold on;

% 10,5 nuDGGK=0.5
eta    = [0.6,0.7,0.8 0.9];
gamma  = [2.4949 0.23368 0.057792 0.023572];
std_g  = [0.8696 0.085267 0.0060463 0.0046137];
shear = [1.7094 0.55278 0.2054 0.085678];
std_s = [0.2428 0.068223 0.012734 0.012291];
subplot(121)
errorbar(eta,gamma,std_g,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.5$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(5,:)); hold on;
subplot(122)
errorbar(eta,shear,std_s,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.5$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(5,:)); hold on;

% Mix len Ricci 2006
subplot(121)
semilogy([0.5  1.0],[10  1e-1],'--','color',[0,0,0]+0.6, 'DisplayName','Mix. Length');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$\Gamma^\infty_{part}$') 
subplot(122)
grate = [0.2194 0.129 0.05084 0.01346];
semilogy([0.6 0.7 0.8 0.9],grate,'--','color',[0,0,0]+0.6, 'DisplayName','Lin. Max Growth');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$V_E$''') 

%% nuDGGK = 0.1
figure; set(gcf, 'Position',  [100, 100, 1200, 350])
% 6,3
eta =   [0.6        0.7     0.8];
gamma = [0.24321    0.085   0.0367];
std_g  = [0.295      0.05    0.0023];
gbmax  = [1.0        0.21    0.0367];
gbmin  = [0.02       0.04    0.0367];
% errorbar(eta,gamma,std_,'-','DisplayName','$P,J=06,03$, $\nu_{DGGK}=0.1$',...
%     'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(6,:)); hold on;
subplot(121)
plot(eta,gamma,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.1$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(6,:)); hold on;
plot(eta,gbmax,'-.','DisplayName',' ','Color',line_colors(6,:)); hold on;
plot(eta,gbmin,'-.','DisplayName',' ','Color',line_colors(6,:)); hold on;

% 10,5
eta =   [0.5        0.6        0.7      0.8       0.9];
gamma = [12.2       0.19761 0.088 0.04253 0.026037];
std_g   = [4.7      0.21328 0.065 0.009 0.00118];
gbmax  = [12.2      0.8 0.25 0.04253 0.026037];
gbmin  = [12.2      0.02 0.03 0.04253 0.026037];
shear = [0.65361 0.46548 0.30645 0.21123 ];
std_s  = [0.21288 0.10233 0.02886 0.0044664 ];
sbmax  = [1 0.7 0.30645 0.21123 ];
sbmin  = [0.4 0.35 0.30645 0.21123 ];
% errorbar(eta,gamma,std_,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.1$',...
%     'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(7,:)); hold on;
subplot(121)
plot(eta,gamma,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.1$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(7,:)); hold on;
plot(eta,gbmax,'-.','DisplayName',' ','Color',line_colors(7,:)); hold on;
plot(eta,gbmin,'-.','DisplayName',' ','Color',line_colors(7,:)); hold on;
subplot(122)
plot(eta(2:end),shear,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.1$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(7,:)); hold on;
plot(eta(2:end),sbmax,'-.','DisplayName',' ','Color',line_colors(7,:)); hold on;
plot(eta(2:end),sbmin,'-.','DisplayName',' ','Color',line_colors(7,:)); hold on;

% Mix len Ricci 2006
subplot(121)
semilogy([0.5  1.0],[10  1e-1],'--','color',[0,0,0]+0.6, 'DisplayName','Mix. Length');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$\Gamma^\infty_{part}$') 
subplot(122)
grate = [0.236 0.174 0.112 0.053];
semilogy([0.6 0.7 0.8 0.9],grate,'--','color',[0,0,0]+0.6, 'DisplayName','Lin. Max Growth');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$V_E$''') 

%% nuSGGK = 1.0
figure
% 6,3 nuDGGK=1.0
eta = [0.5 0.6,0.7,0.8];
gamma = [2.3 0.2215 0.0133 0.0032];
std_g   = [3.1 0.22 0.0019 0.0006];
errorbar(eta,gamma,std_g,'-','DisplayName','$P,J=06,03$, $\nu_{SGGK}=1.0$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(8,:)); hold on;


% 10,5 nuDGGK = 1.0
eta = [0.5 0.6,0.7,0.8];
gamma = [10 0.319 0.009 0.0026];
std_g   = [1.34 0.228 0.001 0.001];
% errorbar(eta,gamma,err,'-.','DisplayName','$P,J=10,5$, $\nu_{SGGK}=1.0$',...
%     'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(1,:)); hold on;

% Mix len Ricci 2006
semilogy([0.5  1.0],[10  1e-1],'--','color',[0,0,0]+0.6, 'DisplayName','Mix. Length');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$\Gamma^\infty_{part}$') 
