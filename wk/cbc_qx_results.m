%% CBC BENCHMARK
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


%% DIMITS

KN       = 2.22;
KT       = [1.00 0.90 0.80 0.70 0.60 0.50]*6.96;
gm42     = [32.2 18.8 10.5 5.89 1.74 0.00];
gm42_err = [05.2 03.8 2.10 1.66 0.53 0.00];
gm84     = [0.00 13.2 7.66 2.86 0.00 0.00];
gm84_err = [0.00 2.79 1.93 0.55 0.00 0.00];
gne      = [0.00 0.00 0.00 0.00 0.00 0.00];
gne_err  = [0.00 0.00 0.00 0.00 0.00 0.00];

figure
errorbar(KT,gm42./KT,gm42_err./KT,'o-', 'LineWidth',1.5); hold on;
errorbar(KT,gm84./KT,gm84_err./KT,'o-', 'LineWidth',1.5);
errorbar(KT, gne./KT, gne_err./KT,'x-k','LineWidth',1.5);

% set(gca, 'YScale', 'log')

legend('GM (4,2)','GM (8,4)','Gene')
xlabel('$\eta=\kappa_T/\kappa_N$'); ylabel('$\chi_i$');
