resdir = '../testcases/ExB_shear_atomic_tests/';
MVIN  = ['cd ',resdir];
curdir = pwd;
MVOUT = ['cd ',curdir];
%------
if strcmp(curdir(1:5),'/home')
    MPIRUN   = 'mpirun';
else
    MPIRUN = '/opt/homebrew/bin/mpirun';
end
%------
EXECDIR  = [curdir(1:end-2),'bin/'];
% EXECNAME = 'gyacomo23_dp_O1';
% EXECNAME = 'gyacomo23_dp';
% EXECNAME = 'gyacomo23_sp';
% EXECNAME = 'gyacomo23_test';
EXECNAME = 'gyacomo23_debug';
% ------
NP = '1'; PARA = '1 1 1';
% NP = '2'; PARA = '1 2 1';
% NP = '6'; PARA = '1 6 1';
%------
% Picture, const. zonal mode in phi in nonlin term
% INNAME = '0';
% Picture, bckg ExB shear
% INNAME = '1';
% NL with padrallel modes initialization no ExB shear
% INNAME = '2';
% NL with parallel modes initialization with ExB shear, no  factor
% INNAME = '3';
% NL with parallel modes initialization with ExB shear, with factor
% INNAME = '4';
% Mcmillan et al. 2019 no correction
% INNAME = '5'; 
% Mcmillan et al. 2019 correction
% INNAME = '6';
% Single mode single step
INNAME = '7';
%-------
RUN   = [MPIRUN,' -np ',NP,' ',EXECDIR,EXECNAME,' ',PARA,' ',INNAME];
system([MVIN,'; ',RUN,'; ',MVOUT]);
%
fast_analysis
