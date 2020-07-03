%Plot growth rate vs kz
default_plots_options % Script to set up default plot variables
% with a linear fit of the log evolution
gammas = zeros(numel(kr),numel(kz));
shifts = zeros(numel(kr),numel(kz));

if K_RICCI
    factor = sqrt(1+MODEL.tau_i);
    fchar  = '\times(1+\tau)^{1/2}$';
else
    factor = 1;
    fchar  = '$';
end

% Linear fit of log(Napj)
x1    = timeNi;
itmin = ceil(0.9 * numel(timeNi)); %Take a subset of the time evolution

for ikz = 1:numel(kz)
    fit = polyfit(x1(itmin:end),log(abs(Nipj(itmin:end,ikr,ikz))),1);
    gammas(ikr,ikz) = fit(1);
    shifts(ikr,ikz) = fit(2);
end

fig = figure;
linename = ['$(P,J)=(',num2str(GRID.pmaxe),',',num2str(GRID.jmaxe),')$'];
plot(factor*kz,gammas(ikr,:),'DisplayName',linename);

TITLE  = [];
TITLE = [TITLE,'$\eta_n=',num2str(1.0/MODEL.eta_n),'$, '];
TITLE = [TITLE,'$\eta_B=',num2str(MODEL.eta_B),'$, '];
TITLE = [TITLE,   '$\nu=',num2str(MODEL.nu),'$, '];
%TITLE = [TITLE,   '$k_z=',num2str(GRID.kz),'$'];

title(TITLE);
grid on
legend('show')
xlabel(['$k_z',fchar])
ylabel('$\gamma L_\perp/c_{s} $')

%% Saving fig
if SAVEFIG
    FIGNAME = 'gamma_k';
    save_figure;
end
