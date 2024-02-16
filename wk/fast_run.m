%% QUICK RUN SCRIPT
% This script runs a simulation directly from the Matlab framework. 
% It is meant to run only small problems in linear for benchmarking and 
% debugging purposes since it makes Matlab "busy".

%% Set up the paths for the necessary Matlab modules
gyacomodir = pwd;
gyacomodir = gyacomodir(1:end-2);
% mpirun     = 'mpirun';
mpirun     = '/opt/homebrew/bin/mpirun'; % for macos
addpath(genpath([gyacomodir,'matlab']))         % Add matlab folder
addpath(genpath([gyacomodir,'matlab/plot']))    % Add plot folder
addpath(genpath([gyacomodir,'matlab/compute'])) % Add compute folder
addpath(genpath([gyacomodir,'matlab/load']))    % Add load folder
addpath(genpath([gyacomodir,'wk/parameters']))  % Add parameters folder

%% Setup run or load an executable
RUN = 1; % To run or just to load
default_plots_options
% EXECNAME = 'gyacomo23_sp_save'; % single precision
% EXECNAME = 'gyacomo23_sp'; % single precision
EXECNAME = 'gyacomo23_dp'; % double precision
% EXECNAME = 'gyacomo23_debug'; % double precision
%% RUN
setup
% system(['rm fort*.90']);
% Run linear simulation
if RUN
    MVIN =['cd ../results/',SIMID,'/',PARAMS,'/;'];
    RUN  =['time ',mpirun,' -np 2 ',gyacomodir,'bin/',EXECNAME,' 1 2 1 0;'];
   % RUN  =['time ',mpirun,' -np 4 ',gyacomodir,'bin/',EXECNAME,' 1 2 2 0;'];
   % RUN  =['time ',mpirun,' -np 6 ',gyacomodir,'bin/',EXECNAME,' 3 2 1 0;'];
     % RUN  =['time ',mpirun,' -np 8 ',gyacomodir,'bin/',EXECNAME,' 2 2 2 0;'];
    % RUN  =['time ',mpirun,' -np 1 ',gyacomodir,'bin/',EXECNAME,' 1 1 1 0;'];
      % RUN = ['./../../../bin/gyacomo23_sp 0;'];
    MVOUT='cd ../../../wk;';
    system([MVIN,RUN,MVOUT]);
end

%% quick load of basic data
filename = [SIMID,'/',PARAMS,'/']; % Create the filename based on SIMID and PARAMS
LOCALDIR = [gyacomodir,'results/',filename,'/']; % Create the local directory path based on gyacomodir, results directory, and filename
FIGDIR   = LOCALDIR; % Set FIGDIR to the same path as LOCALDIR
% Load outputs from jobnummin up to jobnummax
J0 = 0; J1 = 0;
data = {}; % Initialize data as an empty cell array
% load grids, inputs, and time traces
data = compile_results_low_mem(data,LOCALDIR,J0,J1); 

if 0 % Activate or not
%% plot mode evolution and growth rates
% Load phi
[data.PHI, data.Ts3D] = compile_results_3D(LOCALDIR,J0,J1,'phi');
options.NORMALIZED = 0; 
options.TIME   = data.Ts3D;
 % Time window to measure the growth of kx/ky modes
options.KY_TW  = [0.5 1.0]*data.Ts3D(end);
options.KX_TW  = [0.5 1.0]*data.Ts3D(end);
options.NMA    = 1; % Set NMA option to 1
options.NMODES = 999; % Set how much modes we study
options.iz     = 'avg'; % Compressing z
options.ik     = 1; %
options.GOK2   = 0; % plot gamma/k^2
options.fftz.flag = 0; % Set fftz.flag option to 0
options.FIELD = 'phi';
options.SHOWFIG  = 1;
[fig, wkykx, ekykx] = mode_growth_meter(data,options); % Call the function mode_growth_meter with data and options as input arguments, and store the result in fig
end
