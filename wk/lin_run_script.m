%% QUICK RUN SCRIPT
% This script creates a directory in /results and runs a simulation directly
% from the Matlab framework. It is meant to run only small problems in linear
% for benchmarking and debugging purposes since it makes Matlab "busy".

%% Set up the paths for the necessary Matlab modules
gyacomodir = pwd;
gyacomodir = gyacomodir(1:end-2);
mpirun     = 'mpirun';
% mpirun     = '/opt/homebrew/bin/mpirun'; % for macos
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

%% Setup parameters
% run lin_DTT_HM_rho85
% run lin_DTT_HM_rho98
% run lin_DIIID_LM_rho90
% run lin_DIIID_LM_rho95
% run lin_DIIID_LM_rho95_HEL
% run lin_JET_rho97
% run lin_Entropy
% run lin_ITG
% run lin_KBM
% run lin_RHT
% run lin_Ivanov
rho  = 0.95; TRIANG = 'NT'; READPROF = 1; 
% prof_folder = ['parameters/profiles/DIIID_Austin_et_al_2019/',TRIANG,'/'];
% prof_folder = ['parameters/profiles/DIIID_Oak_Nelson/',TRIANG,'/'];
prof_folder = ['parameters/profiles/DIIID_Oak_Nelson_high_density/',TRIANG,'/'];
run lin_DIIID_data
% run lin_STEP_EC_HD_psi71
% run lin_STEP_EC_HD_psi49
if 1
% Plot the profiles
 plot_params_vs_rho(geom,prof_folder,rho,0.5,1.1,Lref,mref,Bref,READPROF);
end
% SIMID = ['rho_scan_DIIID_AUSTIN_2019/3x2x192x96x32/rho',num2str(rho)];
%% Change parameters
% GEOMETRY = 's-alpha';
DELTA =0.0; 
K_Ni = 0; K_Ne = 0;
% DELTA = 0.0; 
% DELTA = 0.2; 
S_DELTA = DELTA/2;
LY   = 2*pi/0.25;
TMAX = 20;
NY   = 2;
DT   = 0.01;
TAU  = 1; NU = 0.05;
% TAU = 1e-3; K_Ti = K_Ti/2/TAU; NU = 3*NU/8/TAU; ADIAB_E = 1; NA = 1;
% MU_X = 1; MU_Y = 1;
%% RUN
setup
system(['cp ../results/',SIMID,'/',PARAMS,'/fort_00.90 ../results/',SIMID,'/.']);
% system(['rm fort*.90']);
% Run linear simulation
if RUN
    MVIN =['cd ../results/',SIMID,'/',PARAMS,'/;'];
    % RUN  =['time ',mpirun,' -np 2 ',gyacomodir,'bin/',EXECNAME,' 1 2 1 0;'];
   RUN  =['time ',mpirun,' -np 4 ',gyacomodir,'bin/',EXECNAME,' 1 4 1 0;'];
   % RUN  =['time ',mpirun,' -np 4 ',gyacomodir,'bin/',EXECNAME,' 2 2 1 0;'];
     % RUN  =['time ',mpirun,' -np 8 ',gyacomodir,'bin/',EXECNAME,' 2 2 2 0;'];
    % RUN  =['time ',mpirun,' -np 1 ',gyacomodir,'bin/',EXECNAME,' 1 1 1 0;'];
      % RUN = ['./../../../bin/gyacomo23_sp 0;'];
    MVOUT='cd ../../../wk;';
    system([MVIN,RUN,MVOUT]);
end

%% Analysis
% load
filename = [SIMID,'/',PARAMS,'/']; % Create the filename based on SIMID and PARAMS
LOCALDIR = [gyacomodir,'results/',filename,'/']; % Create the local directory path based on gyacomodir, results directory, and filename
FIGDIR   = LOCALDIR; % Set FIGDIR to the same path as LOCALDIR
% Load outputs from jobnummin up to jobnummax
J0 = 0; J1 = 0;
data = {}; % Initialize data as an empty cell array
% load grids, inputs, and time traces
data = compile_results_low_mem(data,LOCALDIR,J0,J1); 

if 1 % Activate or not
%% plot mode evolution and growth rates
% Load phi
[data.PHI, data.Ts3D] = compile_results_3D(LOCALDIR,J0,J1,'phi');
options.NORMALIZED = 0; 
options.TIME   = data.Ts3D;
 % Time window to measure the growth of kx/ky modes
options.KY_TW  = [0.25 1.0]*data.Ts3D(end);
options.KX_TW  = [0.25 1.0]*data.Ts3D(end);
options.NMA    = 1; % Set NMA option to 1
options.NMODES = 999; % Set how much modes we study
options.iz     = 'avg'; % Compressing z
options.ik     = 1; %
options.GOK2   = 0; % plot gamma/k^2
options.fftz.flag = 0; % Set fftz.flag option to 0
options.FIELD = 'phi';
options.SHOWFIG = 1;
[fig, wkykx, ekykx] = mode_growth_meter(data,options);
% %%
% kx = (1:data.grids.Nx/2)'*2*pi/data.fort_00.GRID.Lx;
ky = (1:data.grids.Ny/2)'*2*pi/data.fort_00.GRID.Ly;
gkxky = real(wkykx(2:end,1:data.grids.Nx/2))';
gkxky(isnan(gkxky)) =0;
gkxky(isinf(gkxky)) =0;
figure; plot(ky,gkxky(1,:));
% gkxky(gkxky<0)      =0;
% % gkxky = imgaussfilt(gkxky,1);
% %
% wkxky = imag(wkykx(2:end,1:data.grids.Nx/2))';
% wkxky(isnan(wkxky)) =0;
% wkxky(isinf(wkxky)) =0;
% % wkxky(wkxky<0)      =0;
% % wkxky = imgaussfilt(wkxky,1.5);
% %
% figure; 
% subplot(121)
%     contourf(kx,ky,gkxky',10)
%     % clim(0.5*[0 1]); 
%     % colormap(bluewhitered); colorbar;
%     xlim([0.025 1]);
%     xlabel('$k_x\rho_s$'); ylabel('$k_y\rho_s$')
% subplot(122)
%     contourf(kx,ky,wkxky',10)
%     % clim(1*[0 1]); 
%     % colormap(bluewhitered); colorbar 
%     xlim([0.025 1]);
%     xlabel('$k_x\rho_s$'); ylabel('$k_y\rho_s$')
% % save_figure(data,fig,'.png')
end

if (0 && NZ>4)
%% Ballooning plot
[data.PHI, data.Ts3D] = compile_results_3D(LOCALDIR,J0,J1,'phi');
if data.inputs.BETA > 0
[data.PSI, data.Ts3D] = compile_results_3D(LOCALDIR,J0,J1,'psi');
end
options.time_2_plot = [120];
options.kymodes     = [100];
options.normalized  = 1;
options.PLOT_KP     = 0;
% options.field       = 'phi';
options.SHOWFIG  = 1;
[fig, chi, phib, psib, ~] = plot_ballooning(data,options);
end

if 0
%% RH TEST
[data.PHI, data.Ts3D] = compile_results_3D(LOCALDIR,J0,J1,'phi');
ikx = 2; iky = 1; t0 = 1; t1 = data.Ts3D(end);
[~, it0] = min(abs(t0-data.Ts3D));[~, it1] = min(abs(t1-data.Ts3D));
plt = @(x) squeeze(mean(real(x(iky,ikx,:,it0:it1)),3))./squeeze(mean(real(x(iky,ikx,:,it0)),3));
t_ = data.Ts3D(it0:it1);
TH = 1.635*EPS^1.5 + 0.5*EPS^2+0.36*EPS^2.5; theory = 1/(1+Q0^2*TH/EPS^2);
clr_ = lines(20);
figure
plot(t_, -plt(data.PHI)); hold on;
plot(t_,0.5* exp(-t_*NU)+theory,'--k');
plot([t_(1) t_(end)],theory*[1 1],'-k');
plot([t_(1) t_(end)],-mean(plt(data.PHI))*[1 1],'-g');
xlabel('$t$'); ylabel('$\phi_z(t)/\phi_z(0)$')
title(sprintf('$k_x=$%2.2f, $k_y=$%2.2f',data.grids.kx(ikx),data.grids.ky(iky)))
end

if 0
    %% Geometry
    plot_metric(data);
end

if 0
    %% Compiled plot for lin EM analysis
    [data.PHI, data.Ts3D] = compile_results_3D(LOCALDIR,J0,J1,'phi');
    if data.inputs.BETA > 0
    [data.PSI, data.Ts3D] = compile_results_3D(LOCALDIR,J0,J1,'psi');
    end
    options.time_2_plot = [120];
    options.kymodes     = [100];
    options.normalized  = 1;
    options.PLOT_KP     = 0;
    options.SHOWFIG     = 0;
    options.NORMALIZED = 0; 
    options.TIME   = data.Ts3D;
     % Time window to measure the growth of kx/ky modes
    options.KX_TW  = [0.5 1]*data.Ts3D(end);
    options.KY_TW  = [0.5 1]*data.Ts3D(end);
    options.NMA    = 1; % Set NMA option to 1
    options.NMODES = 999; % Set how much modes we study
    options.iz     = 'avg'; % Compressing z
    options.ik     = 1; %
    options.GOK2   = 0; % plot gamma/k^2
    options.fftz.flag = 0; % Set fftz.flag option to 0
    options.FIELD = 'phi';
    options.SHOWFIG  = 0;
    [~, chi, phib, psib, ~] = plot_ballooning(data,options);
    [fig, wkykx, ekykx]     = mode_growth_meter(data,options);
    kx = (1:data.grids.Nx/2)'*2*pi/data.fort_00.GRID.Lx;
    ky = (1:data.grids.Ny/2)'*2*pi/data.fort_00.GRID.Ly;
    [~,~,R,Z] = plot_metric(data,options);
    figure;
    subplot(3,2,1); plot(chi,real(phib),'-b'); hold on; 
                    plot(chi,imag(phib),'-r'); xlabel('$\chi$'); ylabel('$\phi$')
    subplot(3,2,3); plot(chi,real(psib),'-b'); hold on; 
                    plot(chi,imag(psib),'-r'); xlabel('$\chi$'); ylabel('$\psi$')
    subplot(3,2,5); errorbar(ky,squeeze(real(wkykx(2:end,1))),...
                            squeeze(real(ekykx(2:end,1))),...
                            'ok--','MarkerSize',8,'LineWidth',1.5);  hold on;
                    errorbar(ky,squeeze(imag(wkykx(2:end,1))),...
                        squeeze(imag(ekykx(2:end,1))),...
                        '*k--','MarkerSize',8,'LineWidth',1.5);
                    xlabel('$k_y\rho_s$'); ylabel('$\gamma (o),\,\omega (*)$');
    R = R*geom.R0; Z = Z*geom.R0;
    subplot(1,2,2); plot(R,Z,'-k'); xlabel('$R$'); ylabel('$Z$');
    axis equal
    title(data.paramshort);
end