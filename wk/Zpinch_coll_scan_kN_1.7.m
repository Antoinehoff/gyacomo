%%
if 0
figure

Kn = 1.7;

% SUGAMA DK 4,2
% nu_a   = 1e-2*[1.00 2.00 3.00 4.00 5.00 6.00 7.00];
% Gavg_a = 1e-2*[1.00 1.71 2.18 3.11 4.11 5.20 6.08];
% Gstd_a = 1e-2*[1.78 2.67 2.82 3.08 2.33 1.35 1.43];

% errorbar(nu_a, Gavg_a/Kn, Gstd_a/Kn,'DisplayName','Sugama DK (4,2)'); hold on

% SUGAMA GK 4,2
nu_a   = 1e-2*[1.00 2.00 3.00 4.00 5.00 6.00 7.00 8.00 9.00 10.0];
Gavg_a = [2.54e-2 4.66e-2 6.96e-2 8.98e-2 1.06e-1 1.24e-1 1.43e-1 1.52e-1 1.69e-1  1.09e-1];
Gstd_a = [3.04e-2 1.42e-2 1.56e-2 1.23e-2 1.20e-2 1.57e-2 1.63e-2 2.06e-2 2.14e-02 1.78e-2];

errorbar(nu_a, Gavg_a/Kn, Gstd_a/Kn,'DisplayName','Sugama GK (4,2)'); hold on


% FCGK 4,2
nu_a   = 1e-2*[1.00 2.00 3.00 4.00 5.00 6.00 7.00 10.0];
Gavg_a = [8.57e-2 1.45e-1 2.25e-1 2.87e-1 3.48e-1 4.06e-1 4.51e-1 3.65e-1*Kn];
Gstd_a = [2.07e-2 2.61e-2 2.40e-2 3.46e-2 4.30e-2 5.00e-2 5.11e-2  0];

errorbar(nu_a, Gavg_a/Kn, Gstd_a/Kn,'DisplayName','Coulomb (4,2)'); hold on

% LDGK ii 6,3
nu_a   = 1e-2*[1.00 2.00 3.00 4.00 5.00 6.00 7.00 8.00 9.00];
Gavg_a = [3.86e-2 1.82e-2 3.08e-2 5.24e-2 7.08e-2 8.26e-2 5.78e-2 7.16e-2 7.96e-2];
Gstd_a = [3.52e-2 1.87e-2 2.86e-2 2.79e-2 1.72e-2 2.40e-2 2.46e-2 1.01e-2 1.21e-2];

errorbar(nu_a, Gavg_a/Kn, Gstd_a/Kn,'DisplayName','Landau ii (6,3)'); hold on

% Collisionless
plot([0 1], 0.02343*[1 1],'--k','DisplayName','$\nu=0$');


%
xlim([0 0.1]);
legend('show');
xlabel('$\nu R/c_s$'); ylabel('$\Gamma_x^\infty/\kappa_N$');


end
if 0
%%
figure

nu = 0.1;

% FCGK 4,2
kn_a   = [1.60 1.80 2.00 2.20 2.40];
Gavg_a = [1.11e-1 6.86e-1 3.44e-0 1.12e+1 2.87e+1];
Gstd_a = [7.98e-3 1.10e-1 4.03e-1 2.03e+0 7.36e+0];

errorbar(kn_a, Gavg_a./kn_a, Gstd_a./kn_a,'DisplayName','Coulomb (4,2)'); hold on

% % Collisionless
% plot([0 1], 0.02343*[1 1],'--k','DisplayName','$\nu=0$');


%
xlim([1.6 2.5]);
legend('show');
xlabel('$\nu R/c_s$'); ylabel('$\Gamma_x^\infty/\kappa_N$');
end