SIMID = 'MOLI_Comparison'; % Name of the simulations
addpath(genpath('../matlab')) % ... add 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% outputs options
OUTPUTS.nsave_0d = 0;
OUTPUTS.nsave_1d = 0;
OUTPUTS.nsave_2d = 100;
OUTPUTS.nsave_5d = 0;
OUTPUTS.write_Na00    = '.true.';
OUTPUTS.write_moments = '.true.';
OUTPUTS.write_phi     = '.true.';
OUTPUTS.write_doubleprecision = '.true.';
OUTPUTS.resfile0      = '''results''';
%% Grid parameters
GRID.pmaxe = 5;
GRID.jmaxe = 3;
GRID.pmaxi = 5;
GRID.jmaxi = 3;
GRID.nkr   = 1;
GRID.krmin = 0.;
GRID.krmax = 0.;
GRID.nkz   = 1;
GRID.kzmin = 0.1;
GRID.kzmax = 0.1;
%% Model parameters
MODEL.nu      = 0.1;
MODEL.tau_e   = 1.0;
MODEL.tau_i   = 1.0;
MODEL.sigma_e = 0.0233380;
MODEL.sigma_i = 1.0;
MODEL.q_e     =-1.0;
MODEL.q_i     = 1.0;
MODEL.eta_n   = 1.0;
MODEL.eta_T   = 0.0;
MODEL.eta_B   = 0.5;
MODEL.lambdaD = 0.0;
%% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = '''RK4''';
BASIC.nrun                = 100000;
BASIC.dt                  = 0.01;
BASIC.tmax                = 200;
INITIAL.initback_moments  = 0.01;
INITIAL.initnoise_moments = 0.;
INITIAL.iseed             = 42;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write input file
INPUT = write_fort90(OUTPUTS,GRID,MODEL,INITIAL,TIME_INTEGRATION,BASIC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run HeLaZ
nproc = 1;
EXEC  = ' ../bin/helaz ';
RUN   = ['mpirun -np ' num2str(nproc)];
CMD   = [RUN, EXEC, INPUT];
system(CMD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run MOLI
MOLI_time_solver

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis and basic figures
SAVEFIG = 1;
filename = 'results_00.h5';
%% Nipj

moment = 'Ni00';

kr       = h5read(filename,['/data/var2d/' moment '/coordkr']);
kz       = h5read(filename,['/data/var2d/' moment '/coordkz']);
timeNi     = h5read(filename,'/data/var2d/time');
Nipj     = zeros(numel(timeNi),numel(kr),numel(kz));
for it = 1:numel(timeNi)
    tmp          = h5read(filename,['/data/var2d/', moment,'/', num2str(it,'%06d')]);
    Nipj(it,:,:) = tmp.real + 1i * tmp.imaginary; 
end

%% phi
timephi  = h5read(filename,'/data/var2d/time');
kr       = h5read(filename,'/data/var2d/phi/coordkr');
kz       = h5read(filename,'/data/var2d/phi/coordkz');
phiHeLaZ      = zeros(numel(timephi),numel(kr),numel(kz));
for it = 1:numel(timephi)
    tmp         = h5read(filename,['/data/var2d/phi/' num2str(it,'%06d')]);
    phiHeLaZ(it,:,:) = tmp.real + 1i * tmp.imaginary;
end

timephiMOLI = results.timeRK4;
phiMOLI  = zeros(size(timephiMOLI));
for it = 1:numel(timephiMOLI)
    phiMOLI(it) = get_phi(results.NapjRK4(it,:),params,options);
end


%% Error
nsave = OUTPUTS.nsave_2d;
if(numel(phiMOLI(1:nsave:end)) == numel(phiHeLaZ))
    errphi  = abs(phiMOLI(1:nsave:end)-phiHeLaZ)./abs(phiMOLI(1:nsave:end));
    errNipj = abs(results.NapjRK4(1:nsave:end,1)-Nipj)./abs(results.NapjRK4(1:nsave:end,1));
    figure
    plot(timephi,errphi*100,'-','DisplayName','$\epsilon(\phi)$')
    hold on;
    plot(timephi,errNipj*100,'--','DisplayName','$\epsilon(N_i^{00})$')
    title(TITLE);
    xlabel('$t$');
    ylabel('$\epsilon$ [\%]')
    grid on
    legend('show')
else
    figure
    %HeLaZ results
    plot(timephi,abs(phiHeLaZ),'-','DisplayName','HeLaZ RK4')
    title(TITLE);
    hold on
    %MOLI results
    plot(timephiMOLI,abs(phiMOLI),'--','DisplayName','MOLI RK4')
    grid on
    xlabel('$t$')
    ylabel('$|\phi|$')
    legend('show')

    figure
    %HeLaZ results
    x1 = timeNi;
    y1 = abs(Nipj);
    plot(x1,y1,'-','DisplayName','HeLaZ RK4')
    hold on
    %MOLI results
    x2 = results.timeRK4;
    y2 = abs(results.NapjRK4(:,1));
    plot(x2(1:100:end),y2(1:100:end),'--','DisplayName','MOLI RK4');
    title(TITLE);
    grid on
    legend('show')
    xlabel('$t$')
    ylabel(['$|',moment,'|$'])

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%