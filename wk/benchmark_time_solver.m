% Run linear/nonlin simulation on a kr,kz grid and compare with linear MOLI 
clear all; close all;
SIMID = 'non_lin_benchmark_time_solver'; % Name of the simulations
addpath(genpath('../matlab')) % ... add
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% outputs options
OUTPUTS.nsave_0d = 0;
OUTPUTS.nsave_1d = 0;
OUTPUTS.nsave_2d = 20;
OUTPUTS.nsave_5d = 20;
OUTPUTS.write_Na00    = '.false.';
OUTPUTS.write_moments = '.true.';
OUTPUTS.write_phi     = '.true.';
OUTPUTS.write_non_lin = '.false.';
OUTPUTS.write_doubleprecision = '.true.';
OUTPUTS.resfile0      = '''results''';
OUTPUTS.rstfile0      = '''restart''';
%% Grid parameters
GRID.pmaxe = 8;
GRID.jmaxe = 4;
GRID.pmaxi = 8;
GRID.jmaxi = 4;
GRID.nkr   = 8;
GRID.krmin =-2.0;
GRID.krmax = 2.0;
GRID.nkz   = 8;
GRID.kzmin =-2.0;
GRID.kzmax = 2.0;
GRID.Pad   = 2.0;
%% Model parameters
MODEL.CO      = -2;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
MODEL.NON_LIN = '.true.';   % Non linear term
MODEL.nu      = 0.01; % collisionality nu*L_perp/Cs0 (~10^2 bigger than Ricci 2006)
% temperature ratio T_a/T_e
MODEL.tau_e   = 1.0;
MODEL.tau_i   = 1.0;
% mass ratio sqrt(m_a/m_i)
MODEL.sigma_e = 0.0233380;
MODEL.sigma_i = 1.0;
% charge q_a/e
MODEL.q_e     =-1.0;
MODEL.q_i     = 1.0;
% gradients L_perp/L_x
MODEL.eta_n   = 1.0;        % Density
MODEL.eta_T   = 0.0;        % Temperature
MODEL.eta_B   = 0.5;        % Magnetic
% Debye length
MODEL.lambdaD = 0.0;
%% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = '''RK4''';
BASIC.nrun                = 100000;
BASIC.dt                  = 0.01;
BASIC.tmax                = 50.0;
INITIAL.RESTART           = '.true.';
INITIAL.backup_file      = '''restart''';
INITIAL.only_Na00         = '.true.';
INITIAL.initback_moments  = 1.0e-2;
INITIAL.initnoise_moments = 0.0e-5;
INITIAL.iseed             = 42;
INITIAL.selfmat_file = ...
    ['''../iCa/self_Coll_GKE_0_GKI_0_ESELF_1_ISELF_1_Pmaxe_',num2str(GRID.pmaxe),...
    '_Jmaxe_',num2str(GRID.jmaxe),'_Pmaxi_',num2str(GRID.pmaxi),'_Jmaxi_',...
    num2str(GRID.jmaxi),'_pamaxx_10.h5'''];
INITIAL.eimat_file = ...
    ['''../iCa/ei_Coll_GKE_0_GKI_0_ETEST_1_EBACK_1_Pmaxe_',num2str(GRID.pmaxe),...
    '_Jmaxe_',num2str(GRID.jmaxe),'_Pmaxi_',num2str(GRID.pmaxi),'_Jmaxi_',...
    num2str(GRID.jmaxi),'_pamaxx_10_tau_1.0000_mu_0.0233.h5'''];
INITIAL.iemat_file = ...
    ['''../iCa/ie_Coll_GKE_0_GKI_0_ITEST_1_IBACK_1_Pmaxe_',num2str(GRID.pmaxe),...
    '_Jmaxe_',num2str(GRID.jmaxe),'_Pmaxi_',num2str(GRID.pmaxi),'_Jmaxi_',...
    num2str(GRID.jmaxi),'_pamaxx_10_tau_1.0000_mu_0.0233.h5'''];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run HeLaZ
INPUT = write_fort90(OUTPUTS,GRID,MODEL,INITIAL,TIME_INTEGRATION,BASIC);
nproc = 1;
MAKE  = 'cd ..; make; cd wk';
system(MAKE);
EXEC  = ' ../bin/helaz ';
RUN   = ['mpirun -np ' num2str(nproc)];
CMD   = [RUN, EXEC, INPUT];
system(CMD);
% Load results
filename = [OUTPUTS.resfile0(2:end-1),'_00.h5'];
[ Nipj, pgrid_i, jgrid_i, kr, kz, time5d ] = load_5D_data( filename, 'moments_i' );
[ Nepj, pgrid_e, jgrid_e, ~,  ~,  ~ ]      = load_5D_data( filename, 'moments_e' );
[ phiHeLaZ ,~, ~, time2d ]                 = load_2D_data( filename, 'phi' );
if strcmp(OUTPUTS.write_non_lin,'.true.') && strcmp(MODEL.NON_LIN,'.true.')
[ Sepj, ~,  ~, ~,  ~, ~ ]                  = load_5D_data( filename, 'Sepj' );
[ Sipj, ~,  ~, ~,  ~, ~ ]                  = load_5D_data( filename, 'Sipj' );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%^
%% Run MOLI with same input parameters and a given position on k grid
%ikr = ceil(GRID.nkr/2); ikz = ceil(GRID.nkz/2);
ikr = GRID.nkr; ikz = GRID.nkz;
kr_MOLI = kr(ikr); kz_MOLI = kz(ikz);
params.RK4 = 1;
run ../matlab/MOLI_time_solver
if params.RK4; MOLIsolvername = 'RK4';
else;          MOLIsolvername = 'ode15i';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis and basic figures
default_plots_options
SAVEFIG = 1;
default_plots_options
bare = @(p_,j_) (GRID.jmaxe+1)*p_ + j_ + 1;
bari = @(p_,j_) bare(GRID.pmaxe, GRID.jmaxe) + (GRID.jmaxi+1)*p_ + j_ + 1;
%% Plot moments
fig = figure;
    tH = time5d; tM = results.time;
    % Electrons Real%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(331); ic = 1;
        for ip = 0:1
            for ij = 0:1
                yH = squeeze(real(Nepj(ip+1,ij+1,ikr,ikz,:)));
                plot(tH,yH,'-','DisplayName',...
                    ['HeLaZ $N_e^{',num2str(ip),num2str(ij),'}$'],...
                    'color', line_colors(ic,:))
                hold on
                yM = squeeze(real(results.Napj(:,bare(ip,ij))));
                plot(tM,yM,'--','DisplayName',...
                    ['MOLI $N_e^{',num2str(ip),num2str(ij),'}$'],...
                    'color', 0.7*line_colors(ic,:))
                hold on
                ic = ic + 1;
            end
        end
        grid on
        xlabel('$t$'); ylabel(['Re$(N_e^{pj})$'])
    % Electrons Imag%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(334); ic = 1;
        for ip = 0:1
            for ij = 0:1
                yH = squeeze(imag(Nepj(ip+1,ij+1,ikr,ikz,:)));
                plot(tH,yH,'-','DisplayName',...
                    ['HeLaZ $N_e^{',num2str(ip),num2str(ij),'}$'],...
                    'color', line_colors(ic,:))
                hold on
                yM = squeeze(imag(results.Napj(:,bare(ip,ij))));
                plot(tM,yM,'--','DisplayName',...
                    ['MOLI $N_e^{',num2str(ip),num2str(ij),'}$'],...
                    'color', 0.7*line_colors(ic,:))
                hold on
                ic = ic + 1;
            end
        end
        grid on
        xlabel('$t$'); ylabel(['Im$(N_e^{pj})$'])
    % Error on abs(Ne00)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yH = squeeze(abs(Nepj(1,1,ikr,ikz,:))); % HeLaZ
    yM = squeeze(abs(results.Napj(:,bare(0,0))));
    ERR3  = abs(yM(1:OUTPUTS.nsave_5d:numel(yH)*OUTPUTS.nsave_5d)-yH);
    subplot(337);
        yyaxis left
            semilogy(tH,yH,'-','DisplayName','HeLaZ',...
                'color', line_colors(1,:)); hold on;
            semilogy(tM,yM,'--','DisplayName','MOLI',...
                'color', 'k')
            grid on
            xlabel('$t$')
            ylabel('$|N_e^{00}|$')
            legend('show')
        yyaxis right
            semilogy(tH, ERR3,'color', line_colors(2,:));
            set(gca, 'YScale', 'log')
            ylabel('$e(N_i^{00})$')
    % Ions Real%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(332); ic = 1;
        for ip = 0:1
            for ij = 0:1
                yH = squeeze(real(Nipj(ip+1,ij+1,ikr,ikz,:)));
                plot(tH,yH,'-','DisplayName',...
                    ['HeLaZ $N_i^{',num2str(ip),num2str(ij),'}$'],...
                    'color', line_colors(ic,:))
                hold on
                yM = squeeze(real(results.Napj(:,bari(ip,ij))));
                plot(tM,yM,'--','DisplayName',...
                    ['MOLI $N_i^{',num2str(ip),num2str(ij),'}$'],...
                    'color', 0.7*line_colors(ic,:))
                hold on
                ic = ic + 1;
            end
        end
        grid on
        xlabel('$t$'); ylabel(['Re$(N_i^{pj})$'])
    % Ions Imag%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(335); ic = 1;
        for ip = 0:1
            for ij = 0:1
                yH = squeeze(imag(Nipj(ip+1,ij+1,ikr,ikz,:)));
                plot(tH,yH,'-','DisplayName',...
                    ['HeLaZ $N_i^{',num2str(ip),num2str(ij),'}$'],...
                    'color', line_colors(ic,:))
                hold on
                yM = squeeze(imag(results.Napj(:,bari(ip,ij))));
                plot(tM,yM,'--','DisplayName',...
                    ['MOLI $N_i^{',num2str(ip),num2str(ij),'}$'],...
                    'color', 0.7*line_colors(ic,:))
                hold on
                ic = ic + 1;
            end
        end
        grid on
        xlabel('$t$'); ylabel(['Im$(N_i^{pj})$'])
    % Error on abs(Ni00)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yH = squeeze(abs(Nipj(1,1,ikr,ikz,:))); % HeLaZ
    yM = squeeze(abs(results.Napj(:,bari(0,0))));
    ERR3  = abs(yM(1:OUTPUTS.nsave_5d:numel(yH)*OUTPUTS.nsave_5d)-yH);
    subplot(338);
        yyaxis left
            semilogy(tH,yH,'-','DisplayName','HeLaZ',...
                'color', line_colors(1,:)); hold on;
            semilogy(tM,yM,'--','DisplayName','MOLI',...
                'color', 'k')
            grid on
            xlabel('$t$')
            ylabel('$|N_i^{00}|$')
            legend('show')
        yyaxis right
            semilogy(tH, ERR3,'color', line_colors(2,:));
            set(gca, 'YScale', 'log')
            ylabel('$e(N_i^{00})$')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Phi
    timephiMOLI = results.time;
    phiMOLI     = zeros(size(timephiMOLI));
    phiHeLaZ_rz = squeeze(phiHeLaZ(ikr,ikz,:));
    for it = 1:numel(timephiMOLI)
        phiMOLI(it) = get_phi(results.Napj(it,:),params,options);
    end
    subplot(333)
        plot(time2d,real(phiHeLaZ_rz),'-','DisplayName','HeLaZ RK4')
        hold on
        plot(timephiMOLI,real(phiMOLI),'--','DisplayName',['MOLI ',MOLIsolvername])
        grid on
        xlabel('$t$')
        ylabel('Re$(\phi)$')
    %Imag
    subplot(336)
        plot(time2d,imag(phiHeLaZ_rz),'-','DisplayName','HeLaZ RK4')
        hold on
        plot(timephiMOLI,imag(phiMOLI),'--','DisplayName',['MOLI ',MOLIsolvername])
        grid on
        xlabel('$t$')
        ylabel('Im$(\phi)$')
    % phi error
    ERR1 = abs(real(phiMOLI(1:OUTPUTS.nsave_2d:numel(phiHeLaZ_rz)*OUTPUTS.nsave_2d) - phiHeLaZ_rz));
    ERR2 = abs(imag(phiMOLI(1:OUTPUTS.nsave_2d:numel(phiHeLaZ_rz)*OUTPUTS.nsave_2d) - phiHeLaZ_rz));
    subplot(339);
        semilogy(time2d,ERR1,'-','DisplayName','Real')
        hold on
        semilogy(time2d,ERR2,'--','DisplayName','Imag')
        grid on
        xlabel('$t$')
        ylabel('$e(\phi)$')
        legend('show')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finish and save
TITLE  = ['$(k_r,k_z)=(',num2str(kr(ikr)),',',num2str(kz(ikz)),')$, ', TITLE];
suptitle(TITLE);
if strcmp(MODEL.NON_LIN,'.true.'); LINEARITY = 'nl';
else; LINEARITY = 'lin'; end;
if SAVEFIG; FIGNAME = LINEARITY; save_figure; end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%