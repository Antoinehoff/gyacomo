RUN = 1; % To run or just to load
addpath(genpath('../matlab')) % ... add
default_plots_options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set Up parameters
CLUSTER.TIME  = '99:00:00'; % allocation time hh:mm:ss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHYSICAL PARAMETERS
NU      = 0.1;   % Collision frequency
TAU     = 1.0;    % e/i temperature ratio
K_Ne    = 2.0;   % Density gradient drive
K_Ni    = 2.0;   % Density gradient drive
K_Te    = 0.25*K_Ne;   % Temperature '''
K_Ti    = 0.25*K_Ni;   % Temperature '''
K_E     = 0.0;   % Electrostat '''
SIGMA_E = 0.0233380;   % mass ratio sqrt(m_a/m_i) (correct = 0.0233380)
BETA    = 0;
%% GRID PARAMETERS
NX      = 2;     % real space x-gridpoints
NY      = 80;     %     ''     y-gridpoints
LX      = 120;     % Size of the squared frequency domain
LY      = 120;     % Size of the squared frequency domain
NZ      = 1;      % number of perpendicular planes (parallel grid)
SG      = 0;         % Staggered z grids option
%% GEOMETRY
GEOMETRY= 'Z-pinch';
Q0      = 1.0;    % safety factor
SHEAR   = 0.0;    % magnetic shear
EPS     = 0.0;    % inverse aspect ratio
NEXC    = 1;
NPOL    = 1;
COLL_KCUT = 1.8;
%% TIME PARMETERS
TMAX    = 100;  % Maximal time unit
DT      = 1e-2;   % Time step
SPS0D   = 1;      % Sampling per time unit for 2D arrays
SPS2D   = 0;      % Sampling per time unit for 2D arrays
SPS3D   = 1;      % Sampling per time unit for 2D arrays
SPS5D   = 1;    % Sampling per time unit for 5D arrays
SPSCP   = 0;    % Sampling per time unit for checkpoints
JOB2LOAD= -1;
%% OPTIONS
SIMID   = 'Linear_scan';  % Name of the simulation
LINEARITY = 'linear';   % activate non-linearity (is cancelled if KXEQ0 = 1)
KIN_E   = 1;
% Collision operator
% (LB:L.Bernstein, DG:Dougherty, SG:Sugama, LR: Lorentz, LD: Landau)
CO      = 'LR';
GKCO    = 1; % gyrokinetic operator
ABCO    = 1; % INTERSPECIES collisions
INIT_ZF = 0; ZF_AMP = 0.0;
CLOS    = 0;   % Closure model (0: =0 truncation, 1: gyrofluid closure (p+2j<=Pmax))s
NL_CLOS = 0;   % nonlinear closure model (-2:nmax=jmax; -1:nmax=jmax-j; >=0:nmax=NL_CLOS)
KERN    = 0;   % Kernel model (0 : GK)
INIT_OPT= 'mom00';   % Start simulation with a noisy mom00/phi/allmom
%% OUTPUTS
W_DOUBLE = 1;
W_GAMMA  = 1; W_HF     = 1;
W_PHI    = 1; W_NA00   = 1;
W_DENS   = 1; W_TEMP   = 1;
W_NAPJ   = 1; W_SAPJ   = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unused
HD_CO   = 0.0;    % Hyper diffusivity cutoff ratio
kmax    = NX*pi/LX;% Highest fourier mode
MU      = 0.0; % Hyperdiffusivity coefficient
N_HD    = 4;
INIT_BLOB = 0; WIPE_TURB = 0; ACT_ON_MODES = 0;
MU_X    = MU;     % 
MU_Y    = MU;     % 
MU_Z    = 0.0;     %
MU_P    = 0.0;     %
MU_J    = 0.0;     % 
LAMBDAD = 0.0;
NOISE0  = 1.0e-5; % Init noise amplitude
BCKGD0  = 0.0;    % Init background
k_gB   = 1.0;
k_cB   = 1.0;
%% PARAMETER SCANS

if 1
% Parameter scan over PJ
PA = [10];
JA = [5];
Nparam = numel(PA);
% Parameter scan over KN
% PA = [4]; JA = [2];
%     PMAXE = PA(1); PMAXI = PA(1);
%     JMAXE = JA(1); JMAXI = JA(1);
% KNA    = 1.5:0.05:2.5;
% ETA    = 0.25;
% Nparam = numel(KNA);
%
DTA= DT*ones(1,Nparam)./sqrt(JA);
% DTA= DT;
param_name = 'KN';
gamma_Ni00 = zeros(Nparam,numel(ky));
gamma_Nipj = zeros(Nparam,numel(ky));
gamma_phi  = zeros(Nparam,numel(ky));
for i = 1:Nparam
    % Change scan parameter
    PMAXE = PA(i); PMAXI = PA(i);
    JMAXE = JA(i); JMAXI = JA(i);
%     K_N = KNA(i); K_T = ETA*K_N;
    DT = DTA(i);
    setup
    system(['rm fort*.90']);
    % Run linear simulation
    if RUN
        system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 4 ./../../../bin/gyacomo 1 4 1 0; cd ../../../wk'])
% disp([param_name,'=',num2str(K_N)]);
% system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 6 ./../../../bin/helaz3 1 6 0 > out.txt; cd ../../../wk']);
%         system(['cd ../results/',SIMID,'/',PARAMS,'/; mpirun -np 2 ./../../../bin/helaz 1 2 0; cd ../../../wk'])
%         system(['cd ../results/',SIMID,'/',PARAMS,'/; ./../../../bin/helaz 0; cd ../../../wk'])
    end
%     Load and process results
    %%
    filename = ['../results/',SIMID,'/',PARAMS,'/outputs_00.h5'];
    load_results
    for iky = 1:numel(ky)
        tend   = max(Ts3D(abs(Ni00(iky,1,1,:))~=0));
        tstart   = 0.6*tend;
        [~,itstart] = min(abs(Ts3D-tstart));
        [~,itend]   = min(abs(Ts3D-tend));
        trange = itstart:itend;
        % exp fit on moment 00
        X_ = Ts3D(trange); Y_ = squeeze(abs(Ni00(iky,1,1,trange)));
        gamma_Ni00(i,iky) = LinearFit_s(X_,Y_);
        % exp fit on phi
        X_ = Ts3D(trange); Y_ = squeeze(abs(PHI(iky,1,1,trange)));
        gamma_phi (i,iky) = LinearFit_s(X_,Y_);
    end
    gamma_Ni00(i,:) = real(gamma_Ni00(i,:));% .* (gamma_Ni00(i,:)>=0.0));
    gamma_Nipj(i,:) = real(gamma_Nipj(i,:));% .* (gamma_Nipj(i,:)>=0.0));
    if 0
    %% Fit verification
    figure;
    for i = 1:1:numel(ky)
        X_ = Ts3D(:); Y_ = squeeze(abs(Ni00(i,1,1,:)));
        semilogy(X_,Y_,'DisplayName',['k_y=',num2str(ky(i))]); hold on;
    end
end

if 1
%% Plot for PJ scan
SCALE = 1;%sqrt(2);
fig = figure; FIGNAME = 'linear_study';
plt = @(x) x;
% subplot(211)
    for i = 1:Nparam
%         colors = line_colors;
        colors = jet(Nparam);
        clr       = colors(mod(i-1,numel(line_colors(:,1)))+1,:);
        linestyle = line_styles(floor((i-1)/numel(line_colors(:,1)))+1);
        plot(plt(SCALE*ky),plt(gamma_phi(i,1:end)),...
            'Color',clr,...
            'LineStyle',linestyle{1},'Marker','^',...
...%             'DisplayName',['$\kappa_N=',num2str(K_N),'$, $\nu_{',CONAME,'}=',num2str(NU),'$, $P=',num2str(PA(i)),'$, $J=',num2str(JA(i)),'$']);
            'DisplayName',[CONAME,', $P,J=',num2str(PA(i)),',',num2str(JA(i)),'$']);
        hold on;
    end
    grid on; xlabel('$k_y\rho_s^{R}$'); ylabel('$\gamma(\phi)L_\perp/c_s$'); xlim([0.0,max(ky)]);
    title(['$\kappa_N=',num2str(K_Ni),'$, $\nu_{',CONAME,'}=',num2str(NU),'$'])
    legend('show'); %xlim([0.01,10])
saveas(fig,[SIMDIR,'/',PARAMS,'/gamma_vs_',param_name,'_',PARAMS,'.fig']);
saveas(fig,[SIMDIR,'/',PARAMS,'/gamma_vs_',param_name,'_',PARAMS,'.png']);
end
end

if 0
%% Plot for KN scan
fig = figure; FIGNAME = 'linear_study';
[Y,X] = meshgrid(kx,KNA);
pclr = pcolor(Y,X,gamma_phi);set(pclr, 'edgecolor','none'); colorbar;
shading interp
% imagesc(kx,KNA,gamma_phi); 
set(gca,'YDir','normal')
colormap(bluewhitered); xlim([kx(2) kx(end)]);
title(['$\gamma^p$, $\nu_{',CONAME,'}=',num2str(NU),'$'])
xlabel('$k_x$'); ylabel('$\kappa_N$');
end
if 0
%% Space time
    [YT,XT] = meshgrid(Ts3D,kx);
    figure;
%     pclr = surf(XT,YT,squeeze(abs(PHI_ST(1,:,:)))); set(pclr, 'edgecolor','none'); colorbar;
%     pclr = pcolor(XT,YT,squeeze(abs(Ni00_ST(1,:,:)))); set(pclr, 'edgecolor','none'); colorbar;
    semilogy(Ts3D(1:TMAX/SPS3D),squeeze(abs(PHI_ST(1,50:5:100,:))));
end
end