if 0
%%
figure

Kn = 1.7;

% SUGAMA DK 4,2
% nu_a   = 1e-2*[1.00 2.00 3.00 4.00 5.00 6.00 7.00];
% Gavg_a = 1e-2*[1.00 1.71 2.18 3.11 4.11 5.20 6.08];
% Gstd_a = 1e-2*[1.78 2.67 2.82 3.08 2.33 1.35 1.43];

% errorbar(nu_a, Gavg_a/Kn, Gstd_a/Kn,'DisplayName','Sugama DK (4,2)'); hold on

% SUGAMA GK 4,2 200x32
nu_a   = 1e-2*[1.00 2.00 3.00 4.00 5.00 6.00 7.00 8.00 9.00 10.0];
Gavg_a = [2.54e-2 4.66e-2 6.96e-2 8.98e-2 1.06e-1 1.24e-1 1.43e-1 1.52e-1 1.69e-1  1.79e-1];
Gstd_a = [3.04e-2 1.42e-2 1.56e-2 1.23e-2 1.20e-2 1.57e-2 1.63e-2 2.06e-2 2.14e-02 1.78e-2];

errorbar(nu_a, Gavg_a/Kn, Gstd_a/Kn,'DisplayName','Sugama GK (4,2) 200x32'); hold on

% SUGAMA GK 4,2 256x64
nu_a   = 1e-2*[1.00 2.00 3.00 4.00 5.00 6.00 7.00 8.00 9.00 10.0];
Gavg_a = [1.99e-2 3.03e-2 4.44e-2 5.87e-2 7.87e-2 1.07e-1 1.22e-1 1.38e-1 1.63e-1 1.75e-1];
Gstd_a = [1.55e-2 2.22e-2 1.12e-2 1.54e-2 2.20e-2 2.76e-2 2.44e-2 3.20e-2 3.12e-2 3.13e-2];

errorbar(nu_a, Gavg_a/Kn, Gstd_a/Kn,'DisplayName','Sugama GK (4,2) 256x64'); hold on


% SUGAMA GK 6,3 200x32
nu_a   = 1e-2*[1.00 2.00 3.00 4.00 5.00 6.00 7.00 8.00 9.00 10.0];
Gavg_a = [2.35e-2 3.57e-2 5.09e-2 6.97e-2 8.65e-2 1.01e-1 1.16e-1 1.33e-1 1.49e-1 1.62e-1];
Gstd_a = [3.76e-2 3.91e-2 2.96e-2 1.57e-2 1.73e-2 1.94e-2 2.21e-2 2.01e-2 1.93e-2 2.24e-2];

errorbar(nu_a, Gavg_a/Kn, Gstd_a/Kn,'DisplayName','Sugama GK (6,3) 200x32'); hold on


% FCGK 4,2 200x32
nu_a   = 1e-2*[1.00 2.00 3.00 4.00 5.00 6.00 7.00 10.0];
Gavg_a = [8.57e-2 1.45e-1 2.25e-1 2.87e-1 3.48e-1 4.06e-1 4.51e-1 3.65e-1*Kn];
Gstd_a = [2.07e-2 2.61e-2 2.40e-2 3.46e-2 4.30e-2 5.00e-2 5.11e-2  0];

errorbar(nu_a, Gavg_a/Kn, Gstd_a/Kn,'DisplayName','Coulomb (4,2) 200x32'); hold on

% % LDGK ii 6,3 200x32
% nu_a   = 1e-2*[1.00 2.00 3.00 4.00 5.00 6.00 7.00 8.00 9.00];
% Gavg_a = [3.86e-2 1.82e-2 3.08e-2 5.24e-2 7.08e-2 8.26e-2 5.78e-2 7.16e-2 7.96e-2];
% Gstd_a = [3.52e-2 1.87e-2 2.86e-2 2.79e-2 1.72e-2 2.40e-2 2.46e-2 1.01e-2 1.21e-2];
% 
% errorbar(nu_a, Gavg_a/Kn, Gstd_a/Kn,'DisplayName','Landau ii (6,3) 200x32'); hold on

% Collisionless
plot([0 1], 0.02343*[1 1],'--k','DisplayName','$\nu=0$');


%
xlim([0 0.1]);
legend('show');
xlabel('$\nu R/c_s$'); ylabel('$\Gamma_x^\infty/\kappa_N$');
set(gca,'Yscale','log');

end
%% KN SCAN (4,2) 200x64, muHD = 0.1 to 1.0, NL = -1, DT=1e-2 nu = 0.1
if 0
%%
figure

nu = 0.1;


% SGGK 4,2
kn_a   = [1.60 1.80 2.00 2.20 2.40];
Gavg_a = [0.0264    0.7829    3.4742   26.3825   73.3567];
Gstd_a = [0.0039    0.1203    0.3892    4.9250   11.5890];
errorbar(kn_a, Gavg_a./kn_a, Gstd_a./kn_a,'DisplayName','Sugama (4,2)'); hold on

% DGGK 4,2
kn_a   = [1.60 1.80 2.00 2.20 2.40];
Gavg_a = [0.0080    0.0843    1.6752   28.0470   69.5804];
Gstd_a = [0.0010    0.1303    0.3329    6.0057   11.4145];
errorbar(kn_a, Gavg_a./kn_a, Gstd_a./kn_a,'DisplayName','Dougherty (4,2)'); hold on

% LRGK 4,2
kn_a   = [1.60 1.80 2.00 2.20 2.40];
Gavg_a = [0.8244    2.7805    6.4115   21.1632   62.2215];
Gstd_a = [0.0621    0.2755    0.6950    3.6594   15.0836];
errorbar(kn_a, Gavg_a./kn_a, Gstd_a./kn_a,'DisplayName','Lorentz (4,2)'); hold on

% FCGK 4,2
kn_a   = [1.60 1.80 2.00 2.20];
Gavg_a = [0.5071    3.1178   10.1663   29.3928];
Gstd_a = [0.0491    0.3237    1.1579    4.0400];
errorbar(kn_a, Gavg_a./kn_a, Gstd_a./kn_a,'DisplayName','Coulomb (4,2)'); hold on

% LDGKii 4,2
% kn_a   = [1.60 1.80 2.00 2.20 2.40];
% Gavg_a = [1.11e-1 6.86e-1 3.44e-0 1.12e+1 2.87e+1];
% Gstd_a = [7.98e-3 1.10e-1 4.03e-1 2.03e+0 7.36e+0];
% errorbar(kn_a, Gavg_a./kn_a, Gstd_a./kn_a,'DisplayName','Landauii (4,2)'); hold on

% % Collisionless
% plot([0 1], 0.02343*[1 1],'--k','DisplayName','$\nu=0$');


%
title('Transport w.r.t. gradient level at $\nu = 0.1$');
xlim([1.6 2.5]);
legend('show');
xlabel('$L_n/R$'); ylabel('$\Gamma_x^\infty R/L_n$');
end
%% KN SCAN (4,2) 200x64, muHD = 0.1 to 1.0, NL = -1, DT=1e-2 nu = 0.01
if 0
%%
figure
Mksz = 10; Lnwd = 2.0;

nu = 0.01;

pnlty = 1;
% SGGK 4,2
kn_a   = [1.60 1.80 2.00 2.20 2.40 1.70 1.90 2.10 2.30 2.50];
Gavg_a = [0.0104 0.1735 0.6430    3.8566   37.0602 pnlty*0.0433    pnlty*0.5181    pnlty*1.7723   pnlty*14.5713   pnlty*85.4020];
Gstd_a = [0.0045 0.0868 0.3587    1.8440    8.6218 pnlty*0.0214    pnlty*0.2620    pnlty*0.9196   pnlty* 4.1692   pnlty*21.5764];
[~,Idx] = sort(kn_a); 
kn_a = kn_a(Idx); Gavg_a = Gavg_a(Idx); Gstd_a = Gstd_a(Idx);
% errorbar(kn_a, Gavg_a./kn_a, Gstd_a./kn_a,'DisplayName','Sugama (4,2)'); hold on
loglog(kn_a, Gavg_a./kn_a,'s','MarkerSize',Mksz,'LineWidth',Lnwd,'DisplayName','Sugama (4,2)'); hold on

% DGGK 4,2
kn_a   = [1.60 1.80 2.00 2.20 2.40 1.70 1.90 2.10 2.30 2.50];
Gavg_a = [0.0010    0.1245    0.5219    1.7153   38.6814 pnlty*0.0020    pnlty*0.3208    pnlty*1.0034    pnlty*2.6310   pnlty*79.4978];
Gstd_a = [0.0002    0.1154    0.5585    0.9062   12.814  pnlty*0.0005    pnlty*0.1853    pnlty*0.5899    pnlty*1.5473   pnlty*20.7551];
[~,Idx] = sort(kn_a); 
kn_a = kn_a(Idx); Gavg_a = Gavg_a(Idx); Gstd_a = Gstd_a(Idx);
% errorbar(kn_a, Gavg_a./kn_a, Gstd_a./kn_a,'DisplayName','Dougherty (4,2)'); hold on
loglog(kn_a, Gavg_a./kn_a,'v','MarkerSize',Mksz,'LineWidth',Lnwd,'DisplayName','Dougherty (4,2)'); hold on

% LRGK 4,2
kn_a   = [1.60 1.80 2.00 2.20 2.40 1.70 1.90 2.10 2.30 2.50];
Gavg_a = [0.1422    0.5389    1.2577    4.4465    50.4417 pnlty*0.2009    pnlty*0.6764    pnlty*1.8266   pnlty*9.0066  pnlty*83.4325];
Gstd_a = [0.0119    0.1684    0.4304    2.8511    4.7634 pnlty*0.0434    pnlty*0.2677    pnlty*0.9563   pnlty*6.9812  pnlty*27.4134];
[~,Idx] = sort(kn_a); 
kn_a = kn_a(Idx); Gavg_a = Gavg_a(Idx); Gstd_a = Gstd_a(Idx);
% errorbar(kn_a, Gavg_a./kn_a, Gstd_a./kn_a,'DisplayName','Lorentz (4,2)'); hold on
loglog(kn_a, Gavg_a./kn_a,'^','MarkerSize',Mksz,'LineWidth',Lnwd,'DisplayName','Lorentz (4,2)'); hold on

plot(1,1);
% FCGK 4,2
kn_a   = [1.60 1.80 2.00 2.20 2.40 1.70 1.90 2.10 2.30 2.50];
Gavg_a = [0.1914    0.8660    2.2324    5.3947   44.2023 pnlty*0.3803    pnlty*1.3516    pnlty*3.5982   pnlty*16.6212   pnlty*88.5276];
Gstd_a = [0.0097    0.1009    0.4804    1.2491    9.8378 pnlty*0.0202    pnlty*0.3688    pnlty*0.8374    pnlty*3.1974   pnlty*19.0415];
[~,Idx] = sort(kn_a); 
kn_a = kn_a(Idx); Gavg_a = Gavg_a(Idx); Gstd_a = Gstd_a(Idx);
% errorbar(kn_a, Gavg_a./kn_a, Gstd_a./kn_a,'DisplayName','Coulomb (4,2)'); hold on
loglog(kn_a, Gavg_a./kn_a,'d','MarkerSize',Mksz,'LineWidth',Lnwd,'DisplayName','Coulomb (4,2)'); hold on


%
% title('Transport w.r.t. gradient level at $\nu = 0.01$');
xlim([1.55 2.55]);
% legend('show');
xlabel('$\kappa_N$'); ylabel('$\Gamma_x^\infty /\kappa_n$');
grid on
end

if 0 
   %% Mixing length argument
    figure;
   nu = 0.0;
   % collisionless
    kn= [1.6  2.0  2.5 ];
    k = [0.65 0.60 0.50];
    g = [0.19 0.44 0.85];
    y = g.^2./k.^3;
    plot(kn,y*(9.6/y(3)),'--k'); hold on;
    
   nu = 0.01;
   % DG
    kn= [1.6  2.0  2.5 ];
    k = [0.38 0.60 0.66];
    g = [0.13 0.50 0.90];
    y = g.^2./k.^3;
    plot(kn,y*(9.6/y(3)),'--r');  
   % SG
    kn= [1.6  2.0  2.5 ];
    k = [0.41 0.63 0.69];
    g = [0.16 0.51 0.89];
    y = g.^2./k.^3;
    plot(kn,y*(9.6/y(3)),'--b');  
  % LR
    kn= [1.6  2.0  2.5 ];
    k = [0.41 0.57 0.60];
    g = [0.18 0.50 0.87];
    y = g.^2./k.^3;
    plot(kn,y*(9.6/y(3)),'--y');  
  % LD
    kn= [1.6  2.0  2.5 ];
    k = [0.38 0.53 0.57];
    g = [0.15 0.47 0.85];
    y = g.^2./k.^3;
    plot(kn,y*(9.6/y(3)),'--g');  
      
end