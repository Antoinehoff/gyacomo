%% Load results np = 24
np_p_24 = [1   2  3  4];
np_kx_24= [24 12  8  6];
el_ti_np_24= zeros(numel(np_p_24),1);
for i_ = 1:numel(np_p_24)
    npp_ = np_p_24(i_); npkx = np_kx_24(i_);

    %% Load from Marconi
    outfile =['/marconi_scratch/userexternal/ahoffman/HeLaZ/results/Marconi_parallel_scaling_2D/',...
        sprintf('%d_%d',npp_,npkx),...
        '_200x100_L_120_P_12_J_5_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_2e-03/out.txt'];
    BASIC.RESDIR = load_marconi(outfile);
    compile_results
    load_params
    el_ti_np_24(i_) = CPUTIME;
end
%% Load results np = 48
np_p_48 = [1   2  3  4  6];
np_kx_48= [48 24 16 12  8];
el_ti_np_48= zeros(numel(np_p_48),1);
for i_ = 1:numel(np_p_48)
    npp_ = np_p_48(i_); npkx = np_kx_48(i_);

    %% Load from Marconi
    outfile =['/marconi_scratch/userexternal/ahoffman/HeLaZ/results/Marconi_parallel_scaling_2D/',...
        sprintf('%d_%d',npp_,npkx),...
        '_200x100_L_120_P_12_J_5_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_2e-03/out.txt'];
    BASIC.RESDIR = load_marconi(outfile);
    compile_results
    load_params
    el_ti_np_48(i_) = CPUTIME;
end

%% Load results np = 72
np_p_72 = [ 2   3];
np_kx_72= [36  24];
el_ti_np_72= zeros(numel(np_p_72),1);
for i_ = 1:numel(np_p_72)
    npp_ = np_p_72(i_); npkx = np_kx_72(i_);

    %% Load from Marconi
    outfile =['/marconi_scratch/userexternal/ahoffman/HeLaZ/results/Marconi_parallel_scaling_2D/',...
        sprintf('%d_%d',npp_,npkx),...
        '_200x100_L_120_P_12_J_5_eta_0.6_nu_1e-01_DGGK_CLOS_0_mu_2e-03/out.txt'];
    BASIC.RESDIR = load_marconi(outfile);
    compile_results
    load_params
    el_ti_np_72(i_) = CPUTIME;
end

%% Plot
figure
plt = @(x) (x/el_ti_np_24(1)-1)*100;
plot(np_p_24,plt(el_ti_np_24),'o--','DisplayName','Ncpu = 24'); hold on;
plot(np_p_48,plt(el_ti_np_48),'o--','DisplayName','Ncpu = 48 ')
plot(np_p_72,plt(el_ti_np_72),'o--','DisplayName','Ncpu = 72 ')
xlabel('Num. proc. p')
ylabel('Variation from 1 24 [$\%$]')
title('CPU time change from 1D paralel')
legend('show')