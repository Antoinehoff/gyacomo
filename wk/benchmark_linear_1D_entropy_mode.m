addpath(genpath('../matlab')) % ... add
default_plots_options
SIMDIR = '../results/benchmarks/1D_linear_entropy_mode/';
% Compile HeLaZ
cd ..
system('make');
cd wk
% Run HeLaZ in sequential
system(['cd ',SIMDIR,';',' ./../../../bin/helaz 0; cd ../../../wk']);
% Run HeLaZ in parallel (discrepencies can occur at marginal growth rate)
% since the random seed will not be the same)
% system(['cd ',SIMDIR,';',' mpirun -np 6 ./../../../bin/helaz 1 6 0; cd ../../../wk'])
% system(['cd ',SIMDIR,';',' mpirun -np 6 ./../../../bin/helaz 2 3 0; cd ../../../wk'])
%compute growth rate
%%
filename = [SIMDIR,'/outputs_00.h5'];
[ PHI, Ts3D, dt3D] = load_3D_data(filename, 'phi');
[Pe, Je, Pi, Ji, kx, ky, z] = load_grid_data(filename);   
gamma_phi  = kx;
for ikx = 2:size(kx)
    tend   = max(Ts3D(abs(PHI(ikx,1,1,:))~=0));
    tstart   = 0.6*tend;
    [~,itstart] = min(abs(Ts3D-tstart));
    [~,itend]   = min(abs(Ts3D-tend));
    trange = itstart:itend;
    % exp fit on phi
    X_ = Ts3D(trange); Y_ = squeeze(abs(PHI(ikx,1,1,trange)));
    gamma_phi (ikx) = LinearFit_s(X_,Y_);
end

%% Plot
SCALE = 1;%sqrt(2);
openfig([SIMDIR,'/benchmark_data.fig']); hold on;
plot(kx,gamma_phi,'-x','DisplayName','new results');
hold on;
legend('show');

