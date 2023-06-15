d_ = load("rodriguez_fernandez_2020_fig3_sparc_ne_prof.txt");
rho_ne = d_(:,1); ne = d_(:,2);
d_ = load("rodriguez_fernandez_2020_fig3_sparc_T_prof.txt");
rho_Te = d_(:,1); Te = d_(:,2);
figure;
plot(rho_ne,ne,'.k','DisplayName','$n_e$'); hold on;
plot(rho_Te,Te,'or','DisplayName','$T_e$');
%%
[~,i1]    = min(abs(rho_ne-0.95));
[~,i2]    = min(abs(rho_ne-0.98));
[~,in975] = min(abs(rho_ne-0.975));

fit_n = polyfit(rho_ne(i1:i2),ne(i1:i2),1);

rho_ = 0.925:0.025:1.0;

plot(rho_,fit_n(1).*rho_+fit_n(2),'-k','DisplayName',['$-\partial_x n/n_{975}=$',num2str(-fit_n(1)/ne(in975))]);

%%
[~,i1]   = min(abs(rho_Te-0.95));
[~,i2]   = min(abs(rho_Te-0.98));
[~,iT975] = min(abs(rho_Te-0.975));

fit_T = polyfit(rho_Te(i1:i2),Te(i1:i2),1);

rho_ = 0.925:0.025:1.0;

plot(rho_,fit_T(1).*rho_+fit_T(2),'-r','DisplayName',['$-\partial_x T/T_{975}=$',num2str(-fit_T(1)/Te(iT975))]);
xline(rho_Te(iT975));
legend('show')