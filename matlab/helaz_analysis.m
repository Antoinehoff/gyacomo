%% HeLaZ data
filename = 'results_00.h5';
default_plots_options % Script to set up default plot variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data
moment = 'Ni00';

kr       = h5read(filename,['/data/var2d/' moment '/coordkr']);
kz       = h5read(filename,['/data/var2d/' moment '/coordkz']);
timeNi   = h5read(filename,'/data/var2d/time');
Nipj     = zeros(numel(timeNi),numel(kr),numel(kz));
for it = 1:numel(timeNi)
    tmp          = h5read(filename,['/data/var2d/', moment,'/', num2str(it,'%06d')]);
    Nipj(it,:,:) = tmp.real + 1i * tmp.imaginary; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot growth rate vs t
gammas = zeros(numel(kr),numel(kz));
shifts = zeros(numel(kr),numel(kz));
% Linear fit of log(Napj)
x1    = timeNi;
itmin = ceil(0.5 * numel(timeNi)); %Take the second half of the time evolution
for ikr = 1:numel(kr)
    for ikz = 1:numel(kz)
        fit = polyfit(x1(itmin:end),log(abs(Nipj(itmin:end,ikr,ikz))),1);
        gammas(ikr,ikz) = fit(1);
        shifts(ikr,ikz) = fit(2);
    end
end

FIGNAME = 'gamma_t';
fig = figure;
for ikr = 1:numel(kr)
    linename = ['$k_r = ',num2str(kr(ikr)),'$'];
    plot(kz,gammas(ikr,:),'DisplayName',linename);
end
TITLE  = [];
TITLE = [TITLE,'$\eta_n=',num2str(1.0/MODEL.eta_n),'$, '];
TITLE = [TITLE,'$\eta_B=',num2str(MODEL.eta_B),'$, '];
TITLE = [TITLE,   '$\nu=',num2str(MODEL.nu),'$, '];
%TITLE = [TITLE,   '$k_z=',num2str(GRID.kz),'$'];

title(TITLE);
grid on
legend('show')
xlabel('$t$')
ylabel(['$|',moment,'|$'])

%% Saving fig
if SAVEFIG
    FIGDIR = ['../results/', SIMID,'/'];
    if ~exist(FIGDIR, 'dir')
       mkdir(FIGDIR)
    end
    FIGNAME = [FIGNAME,'_Pe_',num2str(GRID.pmaxe),'_Je_',num2str(GRID.jmaxe),...
        '_Pi_',num2str(GRID.pmaxi),'_Ji_',num2str(GRID.jmaxi),...
        '_etan_',num2str(MODEL.eta_n),'_etaB_',num2str(MODEL.eta_B),'_nu_',num2str(MODEL.nu)];
    FIGNAME = [FIGDIR, FIGNAME,'.fig'];
    savefig(fig,FIGNAME);
    disp(['Figure saved @ : ',FIGNAME])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Ni00 evolution
FIGNAME = 'Ni00_t';
fig = figure;
%HeLaZ results
x1 = timeNi;
il = 1;
for ikr = 1:numel(kr)
    ic = 1;
    for ikz = 1:2:numel(kz)
        linename = ['$k_r = ',num2str(kr(ikr)),', k_z = ', num2str(kz(ikz)),'$'];
        y1 = abs(Nipj(:,ikr,ikz));
        semilogy(x1,y1,...
            'DisplayName',linename,'Color', line_colors(ic,:), 'LineStyle', '--')
         hold on
        semilogy(x1(itmin:end),exp(gammas(ikr,ikz)*x1(itmin:end) + shifts(ikr,ikz)),...
            'DisplayName',[linename,' fit'],'Color', line_colors(ic,:), 'LineStyle', '-.')
        ic = ic + 1;
    end
    il = il + 1;
end

TITLE  = [];
TITLE = [TITLE,'$\eta_n=',num2str(1.0/MODEL.eta_n),'$, '];
TITLE = [TITLE,'$\eta_B=',num2str(MODEL.eta_B),'$, '];
TITLE = [TITLE,   '$\nu=',num2str(MODEL.nu),'$, '];
%TITLE = [TITLE,   '$k_z=',num2str(GRID.kz),'$'];

title(TITLE);
grid on
legend('show')
xlabel('$t$')
ylabel(['$|',moment,'|$'])

%% Saving fig
if SAVEFIG
    FIGDIR = ['../results/', SIMID,'/'];
    if ~exist(FIGDIR, 'dir')
       mkdir(FIGDIR)
    end
    FIGNAME = [FIGNAME,'_Pe_',num2str(GRID.pmaxe),'_Je_',num2str(GRID.jmaxe),...
        '_Pi_',num2str(GRID.pmaxi),'_Ji_',num2str(GRID.jmaxi),...
        '_etan_',num2str(MODEL.eta_n),'_etaB_',num2str(MODEL.eta_B),'_nu_',num2str(MODEL.nu)];
    FIGNAME = [FIGDIR, FIGNAME,'.fig'];
    savefig(fig,FIGNAME);
    disp(['Figure saved @ : ',FIGNAME])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%