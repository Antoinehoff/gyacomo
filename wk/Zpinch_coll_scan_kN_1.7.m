figure

Kn = 1.7;

% SUGAMA 4,2
nu_a   = 1e-2*[1.00 2.00 3.00 4.00 5.00 6.00 7.00 10.0];
Gavg_a = 1e-2*[1.00 1.71 2.18 3.11 4.11 5.20 6.08  5.59*Kn];
Gstd_a = 1e-2*[1.78 2.67 2.82 3.08 2.33 1.35 1.43  0.0];

errorbar(nu_a, Gavg_a/Kn, Gstd_a/Kn,'DisplayName','Sugama (4,2)'); hold on

% LANDAU 4,2
nu_a   = 1e-2*[1.00 2.00 3.00 4.00 5.00 6.00 7.00 10.0];
Gavg_a = [8.57e-2 1.45e-1 2.25e-1 2.87e-1 3.48e-1 4.06e-1 4.51e-1 3.65e-1*Kn];
Gstd_a = [2.07e-2 2.61e-2 2.40e-2 3.46e-2 4.30e-2 5.00e-2 5.11e-2  0];

errorbar(nu_a, Gavg_a/Kn, Gstd_a/Kn,'DisplayName','Coulomb (4,2)'); hold on

% LANDAU 6,3
nu_a   = 1e-2*[1.00 2.00 3.00 4.00 5.00 6.00 7.00 8.00 9.00];
Gavg_a = [3.86e-2 1.82e-2 3.08e-2 5.24e-2 7.08e-2 8.26e-2 5.78e-2 7.16e-2 7.96e-2];
Gstd_a = [3.52e-2 1.87e-2 2.86e-2 2.79e-2 1.72e-2 2.40e-2 2.46e-2 1.01e-2 1.21e-2];

errorbar(nu_a, Gavg_a/Kn, Gstd_a/Kn,'DisplayName','Coulomb (6,3)'); hold on

% Collisionless
plot([0 1], 0.02343*[1 1],'--k','DisplayName','$\nu=0$');


%
xlim([min(nu_a) max(nu_a)]);
xlabel('$\nu R/c_s$'); ylabel('$\Gamma_x^\infty/\kappa_N$');