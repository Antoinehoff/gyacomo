%% QUICK RUN SCRIPT
% This script creates a directory in /results and runs a simulation directly
% from the Matlab framework. It is meant to run only small problems in linear
% for benchmarking and debugging purposes since it makes Matlab "busy".

%% Set up the paths for the necessary Matlab modules
wkdir = pwd;
gyacomodir = wkdir(1:end-2);
mpirun     = 'mpirun';
% mpirun     = '/opt/homebrew/bin/mpirun'; % for macos
addpath(genpath([gyacomodir,'matlab']))         % Add matlab folder
addpath(genpath([gyacomodir,'matlab/plot']))    % Add plot folder
addpath(genpath([gyacomodir,'matlab/compute'])) % Add compute folder
addpath(genpath([gyacomodir,'matlab/load']))    % Add load folder
addpath(genpath([gyacomodir,'wk/parameters']))  % Add parameters folder

%% Setup run or load an executable
RUN     = 1; % To run or just to load
RERUN   = 0; % rerun if the  data does not exist
default_plots_options
% EXECNAME = 'gyacomo23_sp'; % single precision
EXECNAME = 'gyacomo23_dp'; % double precision

%% Setup parameters
% run lin_DTT_AB_rho85
% run lin_DTT_AB_rho98
% run lin_JET_rho97
% run lin_Entropy
% run lin_ITG
% run lin_RHT
% rho  = 0.95; TRIANG = 'PT'; READPROF = 1; 
% prof_folder = ['parameters/profiles/DIIID_Austin_et_al_2019/',TRIANG,'/'];
% prof_folder = ['parameters/profiles/DIIID_Oak_Nelson/',TRIANG,'/'];
% prof_folder = ['parameters/profiles/DIIID_Oak_Nelson_high_density/',TRIANG,'/'];
% run lin_DIIID_data
run lin_DIIID_LM_rho95

%% Change parameters
% NU   = 1;
% TAU  = 1;
NY   = 2;
DELTA = 0.0; TRIANG = '';
S_DELTA = DELTA/2;
% EXBRATE = 0;
% S_DELTA = min(2.0,S_DELTA);
% SIGMA_E  = 0.023;
% NEXC = 0;
LX   = 120;
%% Scan parameters
SIMID = [SIMID,TRIANG,'_scan'];
P_a   = [2 4 6 8]; 
J_a   = [1 2 3 4];
ky    = 0.3;
DT    = 0.001; TMAX = 40;
CO    = 'DG';
%% Scan loop
% arrays for the result
g_ky = zeros(numel(P_a),numel(J_a));
g_std= g_ky*0;
w_ky = g_ky*0;
w_std= g_ky*0;
j = 1;
for PMAX = P_a
    i = 1;
    for JMAX = J_a
        LY   = 2*pi/ky;
        DTSAVE0D =  1.0;
        DTSAVE3D =  0.5;
        %% RUN
        setup
        % naming
        filename = [SIMID,'/',PARAMS,'/'];
        LOCALDIR  = [gyacomodir,'results/',filename,'/'];
        % check if data exist to run if no data
        data_ = {};
        try
            data_ = compile_results_low_mem(data_,LOCALDIR,00,00);
            Ntime = numel(data_.Ts0D);
        catch
            data_.outfilenames = [];
        end
        if RUN && (RERUN || isempty(data_.outfilenames) || (Ntime < 10))
            MVIN =['cd ',LOCALDIR,';'];
            % RUNG  =['time ',mpirun,' -np 2 ',gyacomodir,'bin/',EXECNAME,' 1 2 1 0;'];
            RUNG  =['time ',mpirun,' -np 4 ',gyacomodir,'bin/',EXECNAME,' 1 2 2 0;'];
            % RUNG  =['time ',mpirun,' -np 8 ',gyacomodir,'bin/',EXECNAME,' 2 2 2 0;'];
            % RUNG  =['time ',mpirun,' -np 1 ',gyacomodir,'bin/',EXECNAME,' 1 1 1 0;'];
            % RUNG = ['./../../../bin/gyacomo23_sp 0;'];
            MVOUT=['cd ',wkdir,';'];
            system([MVIN,RUNG,MVOUT]);
        end
        data_    = compile_results_low_mem(data_,LOCALDIR,00,00);
        [data_.PHI, data_.Ts3D] = compile_results_3D(LOCALDIR,00,00,'phi');
        if numel(data_.Ts0D)>10
            % Load results after trying to run
            filename = [SIMID,'/',PARAMS,'/'];
            LOCALDIR  = [gyacomodir,'results/',filename,'/'];
    
            data_    = compile_results_low_mem(data_,LOCALDIR,00,00);
            [data_.PHI, data_.Ts3D] = compile_results_3D(LOCALDIR,00,00,'phi');
            options.NORMALIZED = 0; 
            options.TIME   = data_.Ts3D;
             % Time window to measure the growth of kx/ky modes
            options.KY_TW  = [0.7 1.0]*data_.Ts3D(end);
            options.KX_TW  = [0.7 1.0]*data_.Ts3D(end);
            options.NMA    = 1; % Set NMA option to 1
            options.NMODES = 999; % Set how much modes we study
            options.iz     = 'avg'; % Compressing z
            options.ik     = 1; %
            options.GOK2   = 0; % plot gamma/k^2
            options.fftz.flag = 0; % Set fftz.flag option to 0
            options.FIELD = 'phi';
            options.SHOWFIG = 0;
            [fig, wkykx, ekykx] = mode_growth_meter(data_,options);
            % [wkykx,ekykx] = compute_growth_rates(data_.PHI(:,:,:,it1:it2),data_.Ts3D(it1:it2));
            g_ky (i,j) = real(wkykx(2,1));
            g_std(i,j) = real(ekykx(2,1));
            w_ky (i,j) = imag(wkykx(2,1));
            w_std(i,j) = imag(ekykx(2,1));
            [gmax, ikmax] = max(g_ky(i,j));

            msg = sprintf('gmax = %2.2f, kmax = %2.2f',gmax,data_.grids.ky(ikmax)); disp(msg);
        end
        i = i + 1;
    end
    j = j + 1;
end

%% take max growth rate among z coordinate
y_ = g_ky + 1i*w_ky; 
e_ = g_std+ 1i*w_std;

%% Save scan results (gamma)
metadata = {};
if(numel(ky_a)>1 || numel(P_a)>1)
    pmin  = num2str(min(P_a));   pmax = num2str(max(P_a));
    jmin  = num2str(min(J_a));   jmax = num2str(max(J_a));
    kymin = num2str(min(ky_a));  kymax= num2str(max(ky_a));
    filename = [num2str(NX),'x',num2str(NZ),...
                '_P_',pmin,'_',pmax,...
                '_J_',jmin,'_',jmax,...
                '_ky_',ky,...
                '_kN_',num2str(K_Ni),...
                '_',CONAME,'_',num2str(NU),'_be_',num2str(BETA),...,
                '_d_',num2str(DELTA),'.mat'];
    metadata.name   = filename;
    metadata.kymin  = ky;
    metadata.title  = [...
        '$k_y=$',num2str(ky),...
        ', $\nu_{',CONAME,'}=$',num2str(NU),...
        ', $\kappa_T=$',num2str(K_Ti),...
        ', $\kappa_N=$',num2str(K_Ni)];
    metadata.par    = [num2str(NX),'x1x',num2str(NZ)];
    metadata.nscan  = 2;
    metadata.s1name = '$P$';
    metadata.s1     = P_a;
    metadata.s2name = '$J$';
    metadata.s2     = J_a;
    metadata.dname  = '$\gamma c_s/R$';
    metadata.data   = y_;
    metadata.err    = e_;
    save([SIMDIR,filename],'-struct','metadata');
    disp(['saved in ',SIMDIR,filename]);
if 1
    gamma = real(metadata.data); g_err = real(metadata.err);
    omega = imag(metadata.data); w_err = imag(metadata.err);
    gamma = gamma.*(gamma>0.025);
    figure
    colors_ = jet(numel(metadata.s2));
    subplot(121)
    for i = 1:numel(metadata.s2)
        errorbar(metadata.s1,gamma(:,i),0*g_err(:,i),'s-',...
        'LineWidth',2.0,...
        'DisplayName',[metadata.s2name,'=',num2str(metadata.s2(i))],...
        'color',colors_(i,:)); 
        hold on;
    end
    xlabel(metadata.s1name); ylabel(metadata.dname);title(metadata.title);
    xlim([metadata.s1(1) metadata.s1(end)]);
    
    subplot(122)
    for i = 1:numel(metadata.s2)
        errorbar(metadata.s1,omega(:,i),w_err(:,i),'s-',...
        'LineWidth',2.0,...
        'DisplayName',[metadata.s2name,'=',num2str(metadata.s2(i))],...
        'color',colors_(i,:)); 
        hold on;
    end
    xlabel(metadata.s1name); ylabel('$\omega R/c_s$');title(metadata.title);
    xlim([metadata.s1(1) metadata.s1(end)]);
    
    colormap(colors_);
    clb = colorbar;
    clim([1 numel(metadata.s2)+1]);
    clb.Ticks=linspace(metadata.s2(1),metadata.s2(end),numel(metadata.s2));
    clb.Ticks    =1.5:numel(metadata.s2)+1.5;
    clb.TickLabels=metadata.s2;
    clb.Label.String = metadata.s2name;
    clb.Label.Interpreter = 'latex';
    clb.Label.FontSize= 18;
end
    % clear metadata tosave
end
