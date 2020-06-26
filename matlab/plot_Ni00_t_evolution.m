%% Plot the time evolution of the firt ion moment
default_plots_options % Script to set up default plot variables

fig = figure;

LEGEND = [];

x1 = timeNi;

for ikz = 1:2:numel(kz)

    linename = ['$k_r = ',num2str(kr(ikr)),'$, ','$k_z = ',num2str(kz(ikz)),'$'];
    y1 = abs(Nipj(:,ikr,ikz));
    semilogy(x1,y1,'DisplayName',linename)
    LEGEND = [LEGEND, linename];
     hold on

end

for ikz = 1:2:numel(kz)

    semilogy(x1(itmin:end),...
        exp(gammas(ikr,ikz)*x1(itmin:end) + shifts(ikr,ikz)),...
        'Color', 'k', 'LineStyle', '--','HandleVisibility','off')

end

LEGEND = [LEGEND, 'fits'];

TITLE  = [];
TITLE = [TITLE,'$\eta_n=',num2str(MODEL.eta_n),'$, '];
TITLE = [TITLE,'$\eta_B=',num2str(MODEL.eta_B),'$, '];
TITLE = [TITLE,   '$\nu=',num2str(MODEL.nu),'$, '];
TITLE = [TITLE, '$(P,J)=(',num2str(GRID.pmaxe),',',num2str(GRID.jmaxe),')$'];
title(TITLE);
grid on
xlabel('$t$')
ylabel(['$|',moment,'|$'])

%% Saving fig
FIGNAME = 'Ni00_t';
if SAVEFIG
    save_figure;
end