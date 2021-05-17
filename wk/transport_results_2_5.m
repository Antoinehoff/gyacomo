%% nuDGGK = 1.0
figure
% 6,3 nuDGGK=1.0
eta   = [0.5,0.6,0.7,0.8];
gamma = [32.6443 3.6895 0.3744 0.056];
std_  = [6.1529 0.7986 0.0493 0.0134];
errorbar(eta,gamma,std_,'-','DisplayName','$P,J=06,03$, $\nu_{DGGK}=1.0$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(1,:)); hold on;

% 10,5 nuDGGK = 1.0
eta = [0.6,0.7,0.8];
gamma = [3.5546 0.3917 0.0421];
std_   = [0.5846 0.0486 0.0203];
errorbar(eta,gamma,std_,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=1.0$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(2,:)); hold on;

% 12,6 nuDGGK = 1.0
eta = [0.6];
gamma = [4.064];
std_   = [0.7964];
errorbar(eta,gamma,std_,'-','DisplayName','$P,J=12,06$, $\nu_{DGGK}=1.0$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(3,:)); hold on;


% Mix len Ricci 2006
semilogy([0.5  1.0],[10  1e-1],'--','color',[0,0,0]+0.6, 'DisplayName','Mix. Length');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$\Gamma^\infty_{part}$') 

%% nuDGGK = 0.5
figure
% 6,3 nuDGGK=0.5
eta = [0.5,0.6,0.7,0.8];
gamma = [20.511 2.6292 0.2353 0.057];
std_   = [3.67 1.13 0.055 0.004];
errorbar(eta,gamma,std_,'-','DisplayName','$P,J=06,03$, $\nu_{DGGK}=0.5$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(4,:)); hold on;

% 10,5 nuDGGK=0.5
eta    = [0.6,0.7,0.8];
gamma  = [2.4949 0.2578 0.058];
std_   = [0.8696 0.1191  0.011];
errorbar(eta,gamma,std_,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.5$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(5,:)); hold on;

% Mix len Ricci 2006
semilogy([0.5  1.0],[10  1e-1],'--','color',[0,0,0]+0.6, 'DisplayName','Mix. Length');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$\Gamma^\infty_{part}$') 

%% nuDGGK = 0.1
figure
% 6,3
eta =   [0.6        0.7     0.8];
gamma = [0.24321    0.085   0.0367];
std_  = [0.295      0.05    0.0023];
bmax  = [1.0        0.21    0.0367];
bmin  = [0.02       0.04    0.0367];
% errorbar(eta,gamma,std_,'-','DisplayName','$P,J=06,03$, $\nu_{DGGK}=0.1$',...
%     'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(6,:)); hold on;
plot(eta,gamma,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.1$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(6,:)); hold on;
plot(eta,bmax,'-.','DisplayName',' ','Color',line_colors(6,:)); hold on;
plot(eta,bmin,'-.','DisplayName',' ','Color',line_colors(6,:)); hold on;

% 10,5
eta =   [0.5        0.6        0.7      0.8       0.9];
gamma = [12.2       0.2133     0.088    0.04253   0.02435];
std_   = [4.7       0.21328    0.065    0.009     0.00386];
bmax  = [12.2       0.8        0.25     0.04253   0.02435];
bmin  = [12.2       0.02       0.03     0.04253   0.02435];
% errorbar(eta,gamma,std_,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.1$',...
%     'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(7,:)); hold on;
plot(eta,gamma,'-','DisplayName','$P,J=10,05$, $\nu_{DGGK}=0.1$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(7,:)); hold on;
plot(eta,bmax,'-.','DisplayName',' ','Color',line_colors(7,:)); hold on;
plot(eta,bmin,'-.','DisplayName',' ','Color',line_colors(7,:)); hold on;
% Mix len Ricci 2006
semilogy([0.5  1.0],[10  1e-1],'--','color',[0,0,0]+0.6, 'DisplayName','Mix. Length');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$\Gamma^\infty_{part}$') 

%% nuSGGK = 1.0
figure
% 6,3 nuDGGK=1.0
eta = [0.5 0.6,0.7,0.8];
gamma = [2.3 0.2215 0.0133 0.0032];
std_   = [3.1 0.22 0.0019 0.0006];
errorbar(eta,gamma,std_,'-','DisplayName','$P,J=06,03$, $\nu_{SGGK}=1.0$',...
    'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(8,:)); hold on;


% 10,5 nuDGGK = 1.0
eta = [0.5 0.6,0.7,0.8];
gamma = [10 0.319 0.009 0.0026];
std_   = [1.34 0.228 0.001 0.001];
% errorbar(eta,gamma,err,'-.','DisplayName','$P,J=10,5$, $\nu_{SGGK}=1.0$',...
%     'LineWidth',2.0,'MarkerSize',5, 'Color', line_colors(1,:)); hold on;

% Mix len Ricci 2006
semilogy([0.5  1.0],[10  1e-1],'--','color',[0,0,0]+0.6, 'DisplayName','Mix. Length');
set(gca, 'YScale', 'log'); grid on; legend('show')
xlabel('$L_n/R$'); ylabel('$\Gamma^\infty_{part}$') 
