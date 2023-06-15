%% Heat flux convergence for kt=6.96 and 5.3
figure
title('s-$\alpha$ turb. heat flux conv.');
% KT 6.96, nuDGDK = 0.05, 128x64x24, Nexc 5 (cyclic BC)
P   = [2    4    6  8];
Qx  = [47.6 46.1 48.9 51.6];
std_= [4.19 4.98 2.78 6.25];
    errorbar(P,Qx,std_/2,'s-r',...
    'LineWidth',2.0,'MarkerSize',8,...
    'DisplayName','KT 6.9, nuDGDK 0.05'); hold on
xlabel('$P$, $J=P/2$'); ylabel('$Q_x$');
% % KT 6.96, nuDGDK = 0.05, 192x96x24, Nexc 5 (dirichlet BC)
% P   = [4    6    8    10  ];
% Qx  = [44.1 48.9 45.5 00.0];
% std_= [9.63 6.13 13.8 0.00];
%     errorbar(P,Qx,std_/2,'o-r',...
%     'LineWidth',2.0,'MarkerSize',8,...
%     'DisplayName','KT 6.9, nuDGDK 0.05'); hold on
% xlabel('$P$, $J=P/2$'); ylabel('$Q_x$');
% KT 6.96, nuDGDK = 0.05, 192x96x24, Nexc 5 (dirichlet BC)
P   = [4    6    8    10  ];
Qx  = [36.7 44.4 45.5 00.0];
std_= [2.67 4.57 13.8 0.00];
    errorbar(P,Qx,std_/2,'o-r',...
    'LineWidth',2.0,'MarkerSize',8,...
    'DisplayName','KT 6.9, nuDGDK 0.05'); hold on
xlabel('$P$, $J=P/2$'); ylabel('$Q_x$');
% KT 5.3, nuDGDK = 0.05, 128x64x24, Nexc 5
P   = [4     6     8     10    12   ];
Qx  = [0 0 0 0 0];
std_= [0 0 0 0 0];
errorbar(P,Qx,std_/2,'o--r',...
    'LineWidth',2.0,'MarkerSize',8,...
    'DisplayName','KT 5.3, nuDGDK 0.05');
xlabel('$P$, $J=P/2$'); ylabel('$Q_x$');

% GENE RESULTS 6.96 128x64x24, Nexc 5
Nvp = [32    16    8     8];
Qx  = [34.53 37.96 1.948 13.0];
std_= [7.830 6.048 0.629 3.09];
errorbar(Nvp,Qx,std_/2,'o-b',...
    'LineWidth',2.0,'MarkerSize',8,...
    'DisplayName','KT 6.9 GENE');
xlabel('$P=N_{v\parallel}$, $J=P/2=N_\mu$'); ylabel('$Q_x$');

% GENE RESULTS 5.3 128x64x24, Nexc 5
% Comment: the result with nvp = 8 is not trusworthy as GENE does not have
% a linear instability with this resolution.. It is then hard to draw
% any conclusion from it.
Nvp = [32    16    8];
Qx  = [0.284 0.000 0.370];
std_= [0.177 0.000 0.140];
errorbar(Nvp,Qx,std_/2,'o--b',...
    'LineWidth',2.0,'MarkerSize',8,...
    'DisplayName','KT 5.3 GENE');
xlabel('$P=N_{v\parallel}$, $J=P/2=N_\mu$'); ylabel('$Q_x$');

% DIMITS
Dimits_result = 7.67*2*2.5;
plot([0 32], Dimits_result*[1 1],'-.k','DisplayName','Dimits CBC');
legend('show');

%% KT scans 9x5x128x64x24
clrs = lines(10);
figure
kT_ = [6.96      6.3       5.8        5.3         4.8];
Qx_ = [59.21     40.1740   31.3627    16.04       1.0900];
std_= [6.7880    7.5513    9.7815     4.166       0.4995];
errorbar(kT_,Qx_,std_/2,'s--','color',clrs(4,:),...
    'LineWidth',2.0,'MarkerSize',8,...
    'DisplayName','(8,4)');
xlabel('$K_T$'); ylabel('$Q_x$');

%% Add Dimits results on current plot
plot([0 500], Dimits_result*[1 1],'-.k');
%% Add Mandell Miller results on current plot
Mandell_result = 7.67*7.5;
plot([0 500], Mandell_result*[1.3 1.3],'-.k');
plot([0 500], Mandell_result*[1 1],'--k');
plot([0 500], Mandell_result*[0.7 0.7],'-.k');

%% Old results (before the 2j-1 factor in mirror term)
% KT 6.96, nuDGDK = 0.05, 128x64x16, Nexc 1
% P   = [2     4     12];
% Qx  = [50.00 46.00 41.00];
% std_= [6.600 2.300 6.600];
%     errorbar(P,Qx,std_/2,'s-r',...
%     'LineWidth',2.0,'MarkerSize',8,...
%     'DisplayName','KT 6.9, nuDGDK 0.05'); hold on
% xlabel('$P$, $J=P/2$'); ylabel('$Q_x$');
% % KT 6.96, nuDGDK = 0.05, 128x64x24, Nexc 5
% P   = [4     6     8     10   ];
% Qx  = [67.62 67.50 59.21 64.17];
% std_= [15.42 20.32 17.25 16.05];
%     errorbar(P,Qx,std_/2,'o-r',...
%     'LineWidth',2.0,'MarkerSize',8,...
%     'DisplayName','KT 6.9, nuDGDK 0.05'); hold on
% xlabel('$P$, $J=P/2$'); ylabel('$Q_x$');
% % KT 5.3, nuDGDK = 0.05, 128x64x24, Nexc 5
% P   = [4     6     8     10    12   ];
% Qx  = [44.10 21.61 16.04 0.558 0.901];
% std_= [10.61 6.952 4.166 0.025 0.122];
% errorbar(P,Qx,std_/2,'o--r',...
%     'LineWidth',2.0,'MarkerSize',8,...
%     'DisplayName','KT 5.3, nuDGDK 0.05');
% xlabel('$P$, $J=P/2$'); ylabel('$Q_x$');