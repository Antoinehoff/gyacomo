%% Heat flux Qi [R/rhos^2/cs]
kN = 2.22;
%-------------- GM ---------------
%(P,J)=(2,1)
kT_Qi_GM_32 = ...
    [...
     13. 1.5e+2 1.6e+1;...%192x96x16x3x2 kymin=0.02
     11. 5.1e+2 3.5e+2;...%192x96x16x3x2 kymin=0.02
     9.0 9.6e+1 3.0e+1;...%192x96x16x3x2 kymin=0.05
     7.0 5.0e+1 6.6e+0;...%192x96x16x3x2 kymin=0.05
     6.0 3.0e+1 4.8e+0;...%192x96x16x3x2 kymin=0.05
     5.0 1.1e+1 9.4e-1;...%192x96x16x3x2 kymin=0.05
     4.5 9.2e+0 1.6e+0;...%192x96x16x3x2 kymin=0.05
    ];
%(P,J)=(4,2)
kT_Qi_GM_53 = ...
    [...
     13. 2.0e+2 1.2e+1;...%96x64x16x3x2  kymin=0.02 (large box)
%      13. 1.1e+2 2.0e+1;...%128x64x16x5x3 kymin=0.02 (large box)
%      13. 1.3e+2 3.5e+1;...%128x64x16x5x3 kymin=0.05
     11. 1.2e+2 1.6e+1;...%96x64x16x3x2  kymin=0.02 (large box)
%      11. 1.6e+2 1.8e+1;...%128x64x16x5x3 kymin=0.02
%      11. 9.7e+1 2.2e+1;...%128x64x16x5x3 kymin=0.05
     9.0 8.3e+1 2.2e+1;...%128x64x16x5x3 kymin=0.05
%      9.0 7.6e+1 2.3e+1;...%256x128x16x3x2 kymin=0.05 (high res)
     7.0 4.6e+1 2.3e+0;...%128x64x16x5x3 kymin=0.05
     6.0 3.7e+1 6.9e+0;...%128x64x16x5x3 kymin=0.05
     5.3 1.9e+1 2.0e+0;...%128x64x16x5x3 kymin=0.05
     5.0 1.3e+1 3.3e+0;...%128x64x16x5x3 kymin=0.05
     4.5 9.3e+0 1.0e+0;...%128x64x16x5x3 kymin=0.05
    ];
%(P,J)=(12,2) or higher
kT_Qi_GM_122 = ...
    [...
     7.0 4.1e+1 6.6e+0;...%192x96x24x13x7 kymin=0.05
     4.5 9.6e-1 1.5e-1;...%128x64x16x13x2 kymin=0.05
    ];
%-------------- GENE ---------------
kT_Qi_GENE = ...
    [...
     13. 2.7e+2 2.2e+1;...%128x64x16x24x12 kymin=0.02 (large box)
%      13. 2.0e+2 6.6e+1;...%128x64x16x24x12 kymin=0.05
     11. 1.9e+2 1.7e+1;...%128x64x16x24x12 kymin=0.02 (large box)
%      11. 3.3e+2 1.6e+2;...%128x64x16x24x12 kymin=0.05
     9.0 1.1e+2 4.2e+1;...%128x64x16x24x12 kymin=0.05
     7.0 4.1e+1 2.1e+1;...%128x64x16x24x12 kymin=0.05
     5.3 1.1e+1 1.8e+1;...%128x64x16x24x12 kymin=0.05
     4.5 1.9e-1 3.0e-2;...%128x64x16x24x12 kymin=0.05
    ];
%% Heat conductivity Xi [Ln/rhoi^2/cs] computed as Xi = Qi/kT/kN
%init
kT_Xi_GM_32  = kT_Qi_GM_32;
kT_Xi_GM_53  = kT_Qi_GM_53;
kT_Xi_GM_122 = kT_Qi_GM_122;
kT_Xi_GENE   = kT_Qi_GENE;
%scale
for i = 2:3
kT_Xi_GM_32 (:,i) = kT_Qi_GM_32 (:,i)./kT_Qi_GM_32 (:,1)./kN;
kT_Xi_GM_53 (:,i) = kT_Qi_GM_53 (:,i)./kT_Qi_GM_53 (:,1)./kN;
kT_Xi_GM_122(:,i) = kT_Qi_GM_122(:,i)./kT_Qi_GM_122(:,1)./kN;
kT_Xi_GENE  (:,i) = kT_Qi_GENE  (:,i)./kT_Qi_GENE  (:,1)./kN;
end
%% Dimits fig 3 data
KT_DIM      = [4.0 4.5 5.0 6.0 7.0 9.0 12. 14. 16. 18.];
LLNL_GK_DIM = [5.0 0.0; ...
               7.0 2.5;...
               9.0 5.0;...
               12. 8.0;... 
               14. 9.0;...
               16. 9.5;...
               18. 10.];
UCOL_GK_DIM = [5.0 0.5;...
               6.0 1.0;...
               7.0 1.5];
GFL__97_DIM = [4.0 4.0;...
               4.5 5.0;...
               5.0 6.5;...
               7.0 8.0;...
               9.0 11.];
GFL__98_DIM = [5.0 1.5;...
               7.0 5.0;...
               9.0 7.5];
%% Plot
msz = 8.0; lwt = 2.0;
figure;
if 1
xye = kT_Xi_GM_32;
errorbar(xye(:,1), xye(:,2),xye(:,3),'>-','DisplayName','192x96x16x3x2',...
    'MarkerSize',msz,'LineWidth',lwt); hold on
xye = kT_Xi_GM_53;
errorbar(xye(:,1), xye(:,2),xye(:,3),'<-','DisplayName','128x64x16x5x3',...
    'MarkerSize',msz,'LineWidth',lwt); hold on
xye = kT_Xi_GM_122;
errorbar(xye(:,1), xye(:,2),xye(:,3),'^-','DisplayName','128x64x16x13x3',...
    'MarkerSize',msz,'LineWidth',lwt); hold on
xye = kT_Xi_GENE;
errorbar(xye(:,1), xye(:,2),xye(:,3),'+-.k','DisplayName','GENE 128x64x16x24x12',...
    'MarkerSize',msz,'LineWidth',lwt); hold on
end
if 1
   plot(LLNL_GK_DIM(:,1),LLNL_GK_DIM(:,2),'dk--','DisplayName','Dimits GK, LLNL'); hold on
   plot(UCOL_GK_DIM(:,1),UCOL_GK_DIM(:,2),'*k','DisplayName','Dimits PIC, U.COL'); 
   plot(GFL__97_DIM(:,1),GFL__97_DIM(:,2),'+k','DisplayName','Dimits GFL 97'); 
   plot(GFL__98_DIM(:,1),GFL__98_DIM(:,2),'sk','DisplayName','Dimits GFL 98'); 
    
end
plot([6.96 6.96],[0 12],'-.r','DisplayName','CBC');
plot([4.0  4.00],[0 12],'-.b','DisplayName','$\kappa_T^{crit}$');
xlabel('$\kappa_T$'); ylabel('$\chi_i$[$L_n/\rho_i^2 v_{thi}$]');
xlim([0 20]);
ylim([0 12]);
title('Dimits et al. 2000, Fig. 3');
legend('show'); grid on;