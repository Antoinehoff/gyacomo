addpath(genpath('../matlab')) % ... add
default_plots_options

%% KINETIC ELECTRON TEST
% Run HeLaZ
cd ..
system('make');
cd wk
SIMDIR = '../results/benchmarks/shearless_molix_Kin_e/';
system(['cd ',SIMDIR,';',' ./../../../bin/helaz 0; cd ../../../wk'])

% Run molix
% cd ../../molix
% system('make');
% system('./bin/molix > out.txt');
% cd ../HeLaZ/wk
% Compare the results with molix at a given time
time_2_plot = 5.0;
[Y_,X_]=molix_plot_phi([SIMDIR,'molix_phi.h5'],time_2_plot);
filename = [SIMDIR,'/outputs_00.h5'];
[ PHI, Ts3D, dt3D] = load_3D_data(filename, 'phi');
[Pe, Je, Pi, Ji, kx, ky, z] = load_grid_data(filename);

plot_phi_ballooning; hold on
plot(X_/pi,real(Y_),'ob');
plot(X_/pi,imag(Y_),'or');
plot(X_/pi,abs(Y_) ,'ok');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADIABATIC ELECTRON TEST
% Run HeLaZ
cd ..
system('make');
cd wk
SIMDIR = '../results/benchmarks/shearless_molix_Adiab_e/';
system(['cd ',SIMDIR,';',' ./../../../bin/helaz 0; cd ../../../wk'])
% Run molix
% cd ../../molix
% system('make');
% system('./bin/molix > out.txt');
% cd ../HeLaZ/wk
% Compare the results with molix at a given time
time_2_plot = 5.0;
[Y_,X_]=molix_plot_phi([SIMDIR,'molix_phi.h5'],time_2_plot);
[ PHI, Ts3D, dt3D] = load_3D_data([SIMDIR,'outputs_00.h5'], 'phi');
[Pe, Je, Pi, Ji, kx, ky, z] = load_grid_data([SIMDIR,'outputs_00.h5']);

plot_phi_ballooning; hold on
plot(X_/pi,real(Y_),'ob');
plot(X_/pi,imag(Y_),'or');
plot(X_/pi,abs(Y_) ,'ok');
