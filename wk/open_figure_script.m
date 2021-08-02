simname_ = '';
fname_='';
%% Marconi output file
fname_='';
fname_='';
fname_='';
fname_='';
fname_='';
fname_='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/v2.7_P_6_J_3/200x100_L_60_P_6_J_3_eta_0.6_nu_1e-03_SGGK_CLOS_0_mu_1e-03/out.txt';
simname_ = fname_(54:end-8);

%%
% simname_ = '';
% simname_ = '';
% simname_ = '';
% simname_ = '';
% simname_ = '';
% simname_ = '';




%%
figname_ = '/fig/ZF_transport_drphi_';
% figname_ = '/fig/space_time_';
% figname_ = '/fig/phi_shear_phase_space_';

path_    = '../results/';

[~,idx] = max(simname_=='x');

params_  = simname_(idx-3:end);


openfig([path_,simname_,figname_,params_,'.fig']);