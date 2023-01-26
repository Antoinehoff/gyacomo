figure
%% KT 6.96, nuDGDK = 0.05, 128x64x24, Nexc 5
P   = [4     6     8     10   ];
Qx  = [67.62 67.50 59.21 64.17];
std = [15.42 20.32 17.25 16.05];
    errorbar(P,Qx,std/2,'s-',...
    'LineWidth',2.0,...
    'DisplayName','KT 6.96, nuDGDK 0.05'); hold on
xlabel('$P$, $J=P/2$'); ylabel('$Q_x$');
%% KT 5.3, nuDGDK = 0.05, 128x64x24, Nexc 5
P   = [4     6     8     10   ];
Qx  = [44.10 21.61 16.04 0.489];
std = [10.61 6.952 4.166 0.061];
    errorbar(P,Qx,std/2,'s-',...
    'LineWidth',2.0,...
    'DisplayName','KT 5.3, nuDGDK 0.05');
xlabel('$P$, $J=P/2$'); ylabel('$Q_x$');

title('GYAC, turb. heat flux conv.');