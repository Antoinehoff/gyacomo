%% eta = 0.6
%nu Gammainf mu
SGGK_transport = [...
    1.0e+0, 5.6e-1, 0.0e+0;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e+00_SGGK_mu_0e+00/
    5.0e-1, 4.0e-1, 0.0e+0;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_5e-01_SGGK_mu_0e+00/ before wiping ZF
    5.0e-1, 5.2e-1, 0.0e+0;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_5e-01_SGGK_mu_0e+00/ after wiping ZF
    2.5e-1, 3.6e-1, 0.0e+0;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_2e-01_SGGK_mu_0e+00/
    1.0e-1, 2.2e-1, 0.0e+0;... 
    1.0e-1, 1.9e-1, 0.0e+0;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e-01_SGGK_mu_0e+00
    7.5e-2, 1.5e-1, 0.0e+0;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_7e-02_SGGK_mu_0e+00/
    5.0e-2, 1.1e-1, 0.0e+0;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_7e-02_SGGK_mu_0e+00/
    1.0e-2, 4.6e-2, 3.2e-2;...
    ];

DGGK_transport = [...
    1.0e+0, 5.7e+00, 0.0e+0;... % HeLaZ 2.8 P,J=10,5
    5.0e-1, 1.1e+01, 0.0e+0;... % simulation_B/cw_DGGK
    2.5e-1, 6.4e+00, 0.0e+0;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_2e-01_SGGK_mu_0e+00
    1.0e-1, 3.0e+00, 0.0e+0;... % HeLaZ 2.8 P,J=2,1
    1.0e-1, 2.5e+00, 0.0e+0;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e-01_SGGK_mu_0e+00
    7.5e-2, 1.8e+00, 0.0e+0;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_7e-02_SGGK_mu_0e+00
    1.0e-2, 3.3e-01, 0.0e+0;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e-02_DGGK_mu_5e-04
    1.0e-2, 2.2e-01, 1.0e-4;... % HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e-02_DGGK_mu_3e-03
    1.0e-2, 1.4e-01, 3.0e-3;... % HeLaZ 2.8 P,J=6,3
    ];

figure;
semilogy(SGGK_transport(:,1),SGGK_transport(:,2),'.','MarkerSize',30); hold on;
semilogy(DGGK_transport(:,1),DGGK_transport(:,2),'.','MarkerSize',30);
grid on; xlabel('$\nu$'); ylabel('$\Gamma_x$');
legend(['SGGK';'DGGK'])