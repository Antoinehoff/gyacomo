resdir = '/Users/ahoffmann/gyacomo/results/dev/ExB_SF/shearing_picture/';
MVIN  = ['cd ',resdir];
% RUN   = './gyacomo23_dp 0';
% Picture, const. zonal mode in phi in nonlin term
% RUN   = ['/opt/homebrew/bin/mpirun -np 6 ./gyacomo23_dp 1 6 1 0'];
% Picture, bckg ExB shear
% RUN   = '/opt/homebrew/bin/mpirun -np 6 ./gyacomo23_dp 1 6 1 1';
% RUN   = '/opt/homebrew/bin/mpirun -np 1 ./gyacomo23_dp 1 1 1 1';
% Test with parallel modes initialization (should not trigger nonlinear
% coupling)
% RUN   = '/opt/homebrew/bin/mpirun -np 6 ./gyacomo23_dp 1 6 1 3';
RUN   = '/opt/homebrew/bin/mpirun -np 1 ./gyacomo23_dp 1 1 1 3';
% RUN   = '/opt/homebrew/bin/mpirun -np 1 ./gyacomo23_debug 1 1 1 3';
MVOUT = 'cd /Users/ahoffmann/gyacomo/wk';
system([MVIN,'; ',RUN,'; ',MVOUT]);
%%
fast_analysis
