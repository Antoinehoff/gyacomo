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
rho  = 0.95; TRIANG = 'PT'; READPROF = 1; 
% prof_folder = ['parameters/profiles/DIIID_Austin_et_al_2019/',TRIANG,'/'];
% prof_folder = ['parameters/profiles/DIIID_Oak_Nelson/',TRIANG,'/'];
prof_folder = ['parameters/profiles/DIIID_Oak_Nelson_high_density/',TRIANG,'/'];
run lin_DIIID_data

%% Change parameters
% NU   = 1;
% TAU  = 1;
NY   = 2;
EXBRATE = 0;
S_DELTA = min(2.0,S_DELTA);
SIGMA_E  = 0.023;
NEXC = 0;
LX   = 120;
%% Scan parameters
SIMID = [SIMID,'_scan'];
P_a   = [2 4];
% P_a   = 2;
ky_a  = [0.01 0.02 0.05 0.1  0.2  0.5  1.0  2.0  5.0  10.0];
% ky_a  = 4.0;
dt_a  = logspace(-2,-3,numel(ky_a));
CO    = 'DG';
%% Scan loop
% arrays for the result
g_ky = zeros(numel(ky_a),numel(P_a));
g_std= g_ky*0;
w_ky = g_ky*0;
w_std= g_ky*0;
j = 1;
for PMAX = P_a
    JMAX = P/2;
    i = 1;
    for ky = ky_a
        LY   = 2*pi/ky;
        DT   = dt_a(i);%1e-3;%/(1+log(ky/0.05));%min(1e-2,1e-3/ky);
        TMAX = DT*10000;%2;%min(10,1.5/ky);
        DTSAVE0D = 100*DT;
        DTSAVE3D =  10*DT;
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
if(numel(ky_a)>1 || numel(P_a)>1)
    pmin  = num2str(min(P_a));   pmax = num2str(max(P_a));
    kymin = num2str(min(ky_a));  kymax= num2str(max(ky_a));
    filename = [num2str(NX),'x',num2str(NZ),...
                '_ky_',kymin,'_',kymax,...
                '_P_',pmin,'_',pmax,...
                '_kN_',num2str(K_Ni),...
                '_',CONAME,'_',num2str(NU),'_be_',num2str(BETA),'.mat'];
    metadata.name   = filename;
    metadata.kymin  = ky;
    metadata.title  = ['$\nu_{',CONAME,'}=$',num2str(NU),'$\kappa_T=$',num2str(K_Ti),', $\kappa_N=$',num2str(K_Ni)];
    metadata.par    = [num2str(NX),'x1x',num2str(NZ)];
    metadata.nscan  = 2;
    metadata.s2name = '$P$';
    metadata.s2     = P_a;
    metadata.s1name = '$ky$';
    metadata.s1     = ky_a;
    metadata.dname  = '$\gamma c_s/R$';
    metadata.data   = y_;
    metadata.err    = e_;
    save([SIMDIR,filename],'-struct','metadata');
    disp(['saved in ',SIMDIR,filename]);
    clear metadata tosave
end
