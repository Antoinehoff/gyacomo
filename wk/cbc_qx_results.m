cbc      = [0080 0100 0120];
gm42     = [15.4 32.2 43.2];
gm42_err = [2.22 05.2 08.1];
gm84     = [11.3 25.2 44.7];
gm84_err = [2.45 03.4 06.2];
gne      = [8.44 24.1 43.5];
gne_err  = [1.55 04.3 05.9];
kN    = [1.776 2.22 2.664];
kT    = [5.568 6.96 8.352];
figure
errorbar(cbc,gm42,gm42_err,'o-','LineWidth',1.5); hold on;
errorbar(cbc,gm84,gm84_err,'o-','LineWidth',1.5);
errorbar(cbc,gne,gne_err,'x-k','LineWidth',1.5);

% set(gca, 'YScale', 'log')

legend('GM (4,2)','GM (8,4)','Gene')
xlabel('CBC drive [\%]'); ylabel('Radial Heat Flux $Q_x^\infty$');