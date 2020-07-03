clear all; close all;
SIMID = 'benchmark_time_solver'; % Name of the simulations
addpath(genpath('../matlab')) % ... add 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% outputs options
OUTPUTS.nsave_0d = 0;
OUTPUTS.nsave_1d = 0;
OUTPUTS.nsave_2d = 1;
OUTPUTS.nsave_5d = 1;
OUTPUTS.write_Ni00    = '.true.';
OUTPUTS.write_moments = '.true.';
OUTPUTS.write_phi     = '.true.';
OUTPUTS.write_doubleprecision = '.true.';
OUTPUTS.resfile0      = '''results''';
%% Grid parameters
GRID.pmaxe = 6;
GRID.jmaxe = 6;
GRID.pmaxi = 6;
GRID.jmaxi = 6;
GRID.nkr   = 1;
GRID.krmin = 0.;
GRID.krmax = 0.;
GRID.nkz   = 1;
GRID.kzmin = 1.0;
GRID.kzmax = 1.0;
%% Model parameters
MODEL.CO      = -2;  % Collision operator (0 : L.Bernstein, -1 : Full Coulomb, -2 : Dougherty)
MODEL.nu      = 1.0; % collisionality nu*L_perp/Cs0
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
MODEL.eta_n   = 0.0;        % Density
MODEL.eta_T   = 0.0;        % Temperature
MODEL.eta_B   = 0.0;        % Magnetic
% Debye length
MODEL.lambdaD = 0.0;
%% Time integration and intialization parameters
TIME_INTEGRATION.numerical_scheme  = '''RK4''';
BASIC.nrun                = 100000;
BASIC.dt                  = 0.01;
BASIC.tmax                = 5.0;
INITIAL.initback_moments  = 0.01;
INITIAL.initnoise_moments = 0.;
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
%% Run MOLI with same input parameters
params.RK4 = 1;
run ../matlab/MOLI_time_solver
if params.RK4; MOLIsolvername = 'RK4';
else;          MOLIsolvername = 'ode15i';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Analysis and basic figures
default_plots_options
SAVEFIG = 1;
filename = 'results_00.h5';
TITLE  = [];
TITLE = [TITLE,'$\eta_n=',num2str(MODEL.eta_n),'$, '];
TITLE = [TITLE,'$\eta_B=',num2str(MODEL.eta_B),'$, '];
TITLE = [TITLE,'$\eta_T=',num2str(MODEL.eta_T),'$, '];
TITLE = [TITLE,   '$\nu=',num2str(MODEL.nu),'$, '];
TITLE = [TITLE, '$(P,J)=(',num2str(GRID.pmaxe),',',num2str(GRID.jmaxe),')$'];
if     MODEL.CO == -1; CONAME = 'FC';
elseif MODEL.CO == -2; CONAME = 'DC';
elseif MODEL.CO ==  0; CONAME = 'LB'; end;

bare = @(pp,jj) (GRID.jmaxe +1)*pp + jj+1;                      % electron 1D-index
bari = @(pp,jj) bare(GRID.pmaxe, GRID.jmaxe)+(GRID.jmaxi +1)*pp + jj+1; % ion 1D-index

%% Nepj
% 
% moment = 'moments_e';
% 
% kr       = h5read(filename,['/data/var5d/' moment '/coordkr']);
% kz       = h5read(filename,['/data/var5d/' moment '/coordkz']);
% time   = h5read(filename,'/data/var5d/time');
% Nepj     = zeros(GRID.pmaxe+1, GRID.jmaxe+1,numel(kr),numel(kz),numel(time));
% for it = 1:numel(time)
%     tmp          = h5read(filename,['/data/var5d/', moment,'/', num2str(it,'%06d')]);
%     Nepj(:,:,:,:,it) = tmp.real + 1i * tmp.imaginary; 
% end
% 
% fig = figure;
% 
% %HeLaZ results
% x1 = time;
% y1 = x1;
% ic = 1;
% for ikr = 1:GRID.nkr
%     for ikz = 1:GRID.nkz
%         for ip = 0:1
%             for ij = 0:1
%                 for it = 1:numel(time)
%                     y1(it) = abs(Nepj(ip+1,ij+1,ikr,ikz,it));
%                 end
%                 semilogy(x1,y1,'-','DisplayName',...
%                     ['HeLaZ $N_e^{',num2str(ip),num2str(ij),'}$'],...
%                     'color', line_colors(ic,:))
%                 hold on
%                 ic = ic + 1;
%             end
%         end
%     end
% end
% 
% %MOLI results
% x2 = results.time;
% y2 = x2;
% ic = 1;
% for ip = 0:1
%     for ij = 0:1
%         for it = 1:numel(x2)
%             y2(it) = abs(results.Napj(it,bare(ip,ij)));
%         end
%         semilogy(x2,y2,'--','DisplayName',...
%             ['MOLI $N_e^{',num2str(ip),num2str(ij),'}$'],...
%             'color', line_colors(ic,:))
%         hold on
%         ic = ic + 1;
%     end
% end
% 
% title(TITLE);
% grid on
% legend('show')
% xlabel('$t$')
% ylabel(['$|N_e^{pj}|$'])
% if SAVEFIG; FIGNAME = ['Nepj_kz_',num2str(GRID.kzmin)]; save_figure; end;
% 
%% Nipj

moment = 'moments_i';

kr       = h5read(filename,['/data/var5d/' moment '/coordkr']);
kz       = h5read(filename,['/data/var5d/' moment '/coordkz']);
time   = h5read(filename,'/data/var5d/time');
Nipj     = zeros(GRID.pmaxe+1, GRID.jmaxe+1,numel(kr),numel(kz),numel(time));
for it = 1:numel(time)
    tmp          = h5read(filename,['/data/var5d/', moment,'/', num2str(it,'%06d')]);
    Nipj(:,:,:,:,it) = tmp.real + 1i * tmp.imaginary; 
end

fig = figure;

x1 = time;
x2 = results.time;
ic = 1;

% Real part
subplot(321)
for ip = 0:1
    for ij = 0:1
        y1 = squeeze(real(Nipj(ip+1,ij+1,1,1,:)));
        plot(x1,y1,'-','DisplayName',...
            ['HeLaZ $N_i^{',num2str(ip),num2str(ij),'}$'],...
            'color', line_colors(ic,:))
        hold on
        y2 = squeeze(real(results.Napj(:,bari(ip,ij))));
        plot(x2,y2,'--','DisplayName',...
            ['MOLI $N_i^{',num2str(ip),num2str(ij),'}$'],...
            'color', line_colors(ic,:))
        hold on
        ic = ic + 1;
    end
end
grid on
xlabel('$t$')
ylabel(['Re$(N_i^{pj})$'])

% Im part
ic = 1;
subplot(322)
for ip = 0:1
    for ij = 0:1
        y1 = squeeze(imag(Nipj(ip+1,ij+1,1,1,:)));
        plot(x1,y1,'-','DisplayName',...
            ['HeLaZ $N_i^{',num2str(ip),num2str(ij),'}$'],...
            'color', line_colors(ic,:))
        hold on
        y2 = squeeze(imag(results.Napj(:,bari(ip,ij))));
        plot(x2,y2,'--','DisplayName',...
            ['MOLI $N_i^{',num2str(ip),num2str(ij),'}$'],...
            'color', line_colors(ic,:))
        hold on
        ic = ic + 1;
    end
end
grid on
xlabel('$t$')
ylabel(['Im$(N_i^{pj})$'])
%suptitle(TITLE);
%if SAVEFIG; FIGNAME = ['Nipj_kz_',num2str(GRID.kzmin)]; save_figure; end;

%% phi
timephi  = h5read(filename,'/data/var2d/time');
kr       = h5read(filename,'/data/var2d/phi/coordkr');
kz       = h5read(filename,'/data/var2d/phi/coordkz');
phiHeLaZ      = zeros(numel(timephi),numel(kr),numel(kz));
for it = 1:numel(timephi)
    tmp         = h5read(filename,['/data/var2d/phi/' num2str(it,'%06d')]);
    phiHeLaZ(it,:,:) = tmp.real + 1i * tmp.imaginary;
end

timephiMOLI = results.time;
phiMOLI  = zeros(size(timephiMOLI));
for it = 1:numel(timephiMOLI)
    phiMOLI(it) = get_phi(results.Napj(it,:),params,options);
end

%fig = figure;
%Real
subplot(323)
plot(timephi,real(phiHeLaZ),'-','DisplayName','HeLaZ RK4')
hold on
plot(timephiMOLI,real(phiMOLI),'--','DisplayName',['MOLI ',MOLIsolvername])
grid on
xlabel('$t$')
ylabel('Re$(\phi)$')
%Imag
subplot(324)
plot(timephi,imag(phiHeLaZ),'-','DisplayName','HeLaZ RK4')
hold on
plot(timephiMOLI,imag(phiMOLI),'--','DisplayName',['MOLI ',MOLIsolvername])
grid on
xlabel('$t$')
ylabel('Im$(\phi)$')
%if SAVEFIG; FIGNAME = ['phi_kz_',num2str(GRID.kzmin)]; save_figure; end;

%% phi error
timephi  = h5read(filename,'/data/var2d/time');
kr       = h5read(filename,'/data/var2d/phi/coordkr');
kz       = h5read(filename,'/data/var2d/phi/coordkz');
phiHeLaZ      = zeros(numel(timephi),numel(kr),numel(kz));
for it = 1:numel(timephi)
    tmp         = h5read(filename,['/data/var2d/phi/' num2str(it,'%06d')]);
    phiHeLaZ(it,:,:) = tmp.real + 1i * tmp.imaginary;
end

timephiMOLI = results.time;
phiMOLI  = zeros(size(timephiMOLI));
for it = 1:numel(timephiMOLI)
    phiMOLI(it) = get_phi(results.Napj(it,:),params,options);
end

ERR1 = abs(real(phiMOLI - phiHeLaZ));
ERR2 = abs(imag(phiMOLI - phiHeLaZ));

%fig = figure;
subplot(325);
plot(timephiMOLI,ERR1,'-','DisplayName','Real')
hold on
plot(timephiMOLI,ERR2,'--','DisplayName','Imag')
%title(TITLE);
grid on
xlabel('$t$')
ylabel('$|e_\phi|$')
legend('show')
%if SAVEFIG; FIGNAME = ['err_phi_kz_',num2str(GRID.kzmin)]; save_figure; end;

suptitle(TITLE);
if SAVEFIG; FIGNAME = ['kz_',num2str(GRID.kzmin)]; save_figure; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%