resdir = '../testcases/ExB_shear_atomic_tests/';
MVIN  = ['cd ',resdir];
curdir = pwd;
MVOUT = ['cd ',curdir];
%------
% MPIRUN = '/opt/homebrew/bin/mpirun';
MPIRUN   = 'mpirun';
%------
EXECDIR  = '/home/ahoffman/gyacomo/bin/';
EXECNAME = 'gyacomo23_dp_O1';
% EXECNAME = 'gyacomo23_dp';
% EXECNAME = 'gyacomo23_sp';
% EXECNAME = 'gyacomo23_debug';
%------
% NP = '1'; PARA = '1 1 1';
% NP = '2'; PARA = '1 2 1';
NP = '6'; PARA = '1 6 1';
%------
% Picture, const. zonal mode in phi in nonlin term
% INNAME = '0';
% Picture, bckg ExB shear
% INNAME = '1';
% Linear with parallel modes initialization + ExB shear
% INNAME = '2';
% NL with parallel modes initialization no ExB shear
INNAME = '3';
% Ultra reduced NL with parallel modes and 8x8 grid no shear
% INNAME = '4';
% NL with parallel modes initialization + ExB shear
% INNAME = '5';

%-------
RUN   = [MPIRUN,' -np ',NP,' ',EXECDIR,EXECNAME,' ',PARA,' ',INNAME];
system([MVIN,'; ',RUN,'; ',MVOUT]);
%%
fast_analysis
