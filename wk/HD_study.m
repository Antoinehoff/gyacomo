red_ = [0.6350 0.0780 0.1840];
gre_ = [0.4660 0.6740 0.1880];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOUGHERTY
% nu mu
% bursts simulations
b_=[...
    5e-1, 5e+0;...
    1e-1, 1e-1;...
    1e-1, 3e-2;...
    1e-1, 2e-2;...
    1e-1, 1e-2;...
    5e-2, 3e-3;...
    5e-2, 2e-3;...
    1e-2, 1e-2;... % v2.7_P_6_J_3/200x100_L_120_P_6_J_3_eta_0.6_nu_1e-02_DGGK_CLOS_0_mu_1e-02/
    1e-2, 3e-3;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_5e-02_DGGK_mu_3e-03
    1e-2, 5e-4;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e-02_DGGK_mu_5e-04
    1e-2, 1e-4;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e-02_DGGK_mu_3e-03
    1e-3, 2e-2;... % v2.7_P_10_J_5/200x100_L_120_P_10_J_5_eta_0.6_nu_1e-03_DGGK_CLOS_0_mu_2e-02/
    ];
% converged turb plateau simulations
cp_=[...
    1e+0, 1e-2;...
    1e+0, 3e-3;...
    1e+0, 2e-3;...
    5e-1, 2e-3;...
    5e-1, 2e-2;...
    1e-1, 2e-3;...
    1e-1, 1e-3;...
    5e-2, 1e-3;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e-01_DGGK_mu_3e-02
    5e-2, 5e-4;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_5e-02_DGGK_mu_5e-04
    ];
% moving/no turb plateau
dp_=[...
    1e+0, 5e+0;...
    1e+0, 2e+0;...
    1e+0, 1e+0;...
    1e-3, 2e-2;...
    1e-3, 1e-2;...
    1e-3, 5e-3;...
    1e-3, 2.5e-3;...
    1e-3, 1e-3;... %v2.7_P_6_J_3/200x100_L_60_P_6_J_3_eta_0.6_nu_1e-03_DGGK_CLOS_0_mu_1e-03/
    ];
% not sure
rp_=[...
    0,0;...
    ];
figure; set(gcf, 'Position',  [100, 100, 900, 400])
title('Hyperdiffusion study, Dougherty GK')
grid on; xlim([5e-4,5e0]); ylim([5e-5,5e+0]);
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
xlabel('$\nu$'); ylabel('$\mu_{HD}$'); hold on;

% Trajectory of simulations
if 0
% HD_study/200x100_L_200_P_2_J_1_eta_0.6_nu_1e+00_DGGK_CLOS_0_mu_1e-02/
mu_ = [0  1 2 5 5];
nu_ = [1  1 1 1 0.5];
plot(nu_,mu_,'x--','DisplayName','N=200, L=050, P,J=2,1');

% HD_study/300x150_L_200_P_2_J_1_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_1e-02/
mu_ = [1e-2 1e-3];
nu_ = [1e-1 1e-1];
plot(nu_,mu_,'x--','DisplayName','N=300, L=100, P,J=2,1');

% HD_study/300x150_L_200_P_2_J_1_eta_0.6_nu_5e-01_DGGK_CLOS_0_mu_1e-02/
mu_ = [2e-3 2e-2];
nu_ = [5e-1 5e-1];
plot(nu_,mu_,'x--','DisplayName','N=300, L=100, P,J=2,1');

% HD_study/100x50_L_50_P_2_J_1_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_1e-02/
mu_ = [1e-2 5e-3 2e-3];
nu_ = [1e-1 1e-1 1e-1];
plot(nu_,mu_,'x-','DisplayName','N=100, L=050, P,J=2,1');

% HD_study/150x75_L_100_P_2_J_1_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_3e-02/
mu_ = [3e-2 1e-3    0];
nu_ = [1e-1 1e-1 1e-1];
plot(nu_,mu_,'x--','DisplayName','N=150, L=100, P,J=2,1');

% HD_study/150x75_L_100_P_2_J_1_eta_0.6_nu_1e-02_DGGK_CLOS_0_mu_3e-02/
mu_ = [3e-3 1e-3 1e-4];
nu_ = [1e-2 1e-2 1e-2];
plot(nu_,mu_,'x--','DisplayName','N=150, L=100, P,J=2,1');

% HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e-02_DGGK_CLOS_0_mu_3e-02/
mu_ = [3e-3 5e-4 1e-4];
nu_ = [1e-2 1e-2 1e-2];
plot(nu_,mu_,'x--','DisplayName','N=150, L=100, P,J=4,2');

% HD_study/150x75_L_100_P_2_J_1_eta_0.6_nu_5e-02_DGGK_CLOS_0_mu_3e-02/
mu_ = [3e-3 1e-3];
nu_ = [5e-2 5e-2];
plot(nu_,mu_,'x--','DisplayName','N=150, L=100, P,J=2,1');

% v2.7_P_6_J_3/200x100_L_120_P_6_J_3_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_2e-02/
mu_ = [2e-2];
nu_ = [1e-1];
% plot(nu_,mu_,'x--','DisplayName','N=200, L=120, P,J=6,3');
end
scatter( b_(:,1), b_(:,2),'o',...
    'MarkerFaceColor',red_,'MarkerEdgeColor',[0 0 0],'SizeData',50,...
    'DisplayName','Bursts'); 
scatter(cp_(:,1),cp_(:,2),'s',...
    'MarkerFaceColor',gre_,'MarkerEdgeColor',[0 0 0],'SizeData',50,...
    'DisplayName','Converged Plateau'); 
scatter(dp_(:,1),dp_(:,2),'d',...
    'MarkerFaceColor',gre_,'MarkerEdgeColor',[0 0 0],'SizeData',50,...
    'DisplayName','Moving Plateau'); 
scatter(rp_(:,1),rp_(:,2),'h',...
    'MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerEdgeColor',[0 0 0],'SizeData',50,...
    'DisplayName','not sure'); 
plot(0,0,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0], 'DisplayName','$\mu=0$');
legend('show','Location','NorthWest')
scatter(1,5e-5,80,'v','MarkerFaceColor',gre_,'MarkerEdgeColor',[0 0 0]);
% HD_study/150x75_L_100_P_2_J_1_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_3e-02/
scatter(0.1,5e-5,80,'v','MarkerFaceColor',gre_,'MarkerEdgeColor',[0 0 0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUGAMA
% nu mu
% bursts simulations
b_=[...
    1e+0, 1.0e-2;... % v2.7_P_6_J_3/200x100_L_120_P_6_J_3_eta_0.6_nu_1e+00_SGGK_CLOS_0_mu_1e-02/
    1e+0, 3.0e-3;... % v2.7_P_6_J_3/200x100_L_120_P_6_J_3_eta_0.6_nu_1e+00_SGGK_CLOS_0_mu_1e-02/
    5e-1, 3.0e-2;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_5e-01_SGGK_mu_3e-02
    5e-1, 1.6e-2;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_5e-01_SGGK_mu_3e-02
    5e-1, 0;...      % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_5e-01_SGGK_mu_0e+00/out.txt
    1e-1, 2.0e-2;... % v2.7_P_6_J_3/200x100_L_120_P_6_J_3_eta_0.6_nu_1e-01_SGGK_CLOS_0_mu_2e-02/
    1e-1, 1.6e-2;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e-01_GGK_mu_3e-02/
    1e-2, 5.0e-3;... % v2.7_P_6_J_3/200x100_L_120_P_6_J_3_eta_0.6_nu_1e-02_SGGK_CLOS_0_mu_1e-02/
    1e-2, 3.0e-3;... % v2.7_P_6_J_3/200x100_L_120_P_6_J_3_eta_0.6_nu_1e-02_SGGK_CLOS_0_mu_1e-02/
    ];
% converged turb plateau simulations
cp_=[...
    1e-1, 0;... % v2.7_P_2_J_1/200x100_L_120_P_2_J_1_eta_0.6_nu_1e-01_SGGK_CLOS_0_mu_0e+00/
    ];
% moving/no turb plateau
dp_=[...
    1e-1, 0;...      % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e-01_SGGK_mu_0e+00/out.txt
    1e-2, 1.0e-2;... % v2.7_P_6_J_3/200x100_L_120_P_6_J_3_eta_0.6_nu_1e-02_SGGK_CLOS_0_mu_1e-02/
    1e-2, 5.0e-3;... % v2.7_P_6_J_3/200x100_L_120_P_6_J_3_eta_0.6_nu_1e-02_SGGK_CLOS_0_mu_1e-02/
    1e-2, 3.0e-3;... % v2.7_P_6_J_3/200x100_L_120_P_6_J_3_eta_0.6_nu_1e-02_SGGK_CLOS_0_mu_1e-02/
    1e-3, 1.0e-3;... % v2.7_P_6_J_3/200x100_L_60_P_6_J_3_eta_0.6_nu_1e-03_SGGK_CLOS_0_mu_1e-03/
    ];
% not sure
rp_=[...
    1e-1, 3.0e-2;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e-01_SGGK_mu_3e-02
    5e-2, 1.0e-3;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e-01_DGGK_mu_3e-02
    ];
figure; set(gcf, 'Position',  [100, 100, 900, 400])
title('Hyperdiffusion study, Sugama GK')
grid on; xlim([5e-4,5e0]); ylim([5e-5,5e+0]);
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
xlabel('$\nu$'); ylabel('$\mu_{HD}$'); hold on;

% Trajectory of simulations
if 0
% v2.7_P_2_J_1/200x100_L_120_P_2_J_1_eta_0.6_nu_1e-01_SGGK_CLOS_0_mu_0e+00/
mu_ = [0.0];
nu_ = [0.1];
plot(nu_,mu_,'x--','DisplayName','N=200, L=050, P,J=2,1');

% Trajectory of simulations
% v2.7_P_6_J_3/200x100_L_120_P_6_J_3_eta_0.6_nu_1e-01_SGGK_CLOS_0_mu_2e-02/
mu_ = [0.02];
nu_ = [0.1];
plot(nu_,mu_,'x--','DisplayName','N=200, L=050, P,J=2,1');

mu_ = [0.0];
nu_ = [0.1];
plot(nu_,mu_,'x--','DisplayName','N=200, L=050, P,J=2,1');
end

scatter( b_(:,1), b_(:,2),'o',...
    'MarkerFaceColor',red_,'MarkerEdgeColor',[0 0 0],'SizeData',80,...
    'DisplayName','Bursts'); 
% scatter(cp_(:,1),cp_(:,2),'s',...
%     'MarkerFaceColor',gre_,'MarkerEdgeColor',[0 0 0],'SizeData',60,...
%     'DisplayName','Converged Plateau'); 
scatter(dp_(:,1),dp_(:,2),'d',...
    'MarkerFaceColor',gre_,'MarkerEdgeColor',[0 0 0],'SizeData',60,...
    'DisplayName','Plateau'); 
% scatter(rp_(:,1),rp_(:,2),'h',...
%     'MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerEdgeColor',[0 0 0],'SizeData',60,...
%     'DisplayName','not sure'); 
scatter(0,0,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
    'DisplayName','$\mu=0$');
legend('show','Location','NorthWest')

scatter(0.075,5e-5,80,'v','MarkerFaceColor',gre_,'MarkerEdgeColor',[0 0 0],...
    'DisplayName','$\mu=0$');
scatter(0.1,5e-5,80,'v','MarkerFaceColor',gre_,'MarkerEdgeColor',[0 0 0],...
    'DisplayName','$\mu=0$');
scatter(0.25,5e-5,80,'v','MarkerFaceColor',red_,'MarkerEdgeColor',[0 0 0],...
    'DisplayName','$\mu=0$');
scatter(0.5,5e-5,80,'v','MarkerFaceColor',red_,'MarkerEdgeColor',[0 0 0],...
    'DisplayName','$\mu=0$');
scatter(1.0,5e-5,80,'v','MarkerFaceColor',red_,'MarkerEdgeColor',[0 0 0],...
    'DisplayName','$\mu=0$');
