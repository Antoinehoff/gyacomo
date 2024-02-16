gyacomodir = pwd; gyacomodir = gyacomodir(1:end-2); % get code directory
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add
addpath(genpath([gyacomodir,'wk/parameters'])) % ... add

Nsim = 48;
gr_  = zeros(61,Nsim);
i = 1;

tau_a = logspace(0,-3.5,Nsim);
for i = 1:Nsim
    run lin_Ivanov;
    TAU  = tau_a(i);
    NU   = 0.1*3/2/TAU/4;                 % Collision frequency
    K_Ti = 0.35*2/TAU;     % ion Temperature
    DT   = 5e-3;
    PMAX = 2; JMAX = 1;
    % PMAX = 2; JMAX = 0;
    % PMAX = 0; JMAX = 1;
    run fast_run;
    % growth rates
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
    options.SHOWFIG  = 0;
    [fig, wkykx, ekykx] = mode_growth_meter(data,options); % Call the function mode_growth_meter with data and options as input arguments, and store the result in fig
    gr_(:,i) = wkykx(:,1);
end

%
figure
contourf(data.grids.ky,tau_a,real(gr_)',20)
set(gca,"YScale","log")