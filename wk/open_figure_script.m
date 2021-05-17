simname_ = '';

%%
% simname_ = 'v2.5_P_10_J_5/200x100_L_120_P_10_J_5_eta_0.7_nu_1e+00_SGGK_CLOS_0_mu_2e-02';
% simname_ = 'v2.5_P_6_J_3/200x100_L_120_P_6_J_3_eta_0.8_nu_1e+00_SGGK_CLOS_0_mu_2e-02';

% simname_ = 'v2.5_P_2_J_1/200x100_L_120_P_2_J_1_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_2e-02';
% simname_ = 'v2.5_P_10_J_5/200x100_L_120_P_10_J_5_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_2e-02';
simname_ = 'v2.5_P_6_J_3/200x100_L_120_P_6_J_3_eta_0.6_nu_5e-01_DGGK_CLOS_0_mu_2e-02';


%%
% figname_ = '/fig/ZF_transport_drphi_';
% figname_ = '/fig/space_time_';
figname_ = '/fig/phi_shear_phase_space_';

path_    = '../results/';

[~,idx] = max(simname_=='x');

params_  = simname_(idx-3:end);


openfig([path_,simname_,figname_,params_,'.fig']);