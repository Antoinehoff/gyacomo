
[t_PT, Gi_PT, Qi_PT, Ge_PT, Qe_PT] = read_flux_out_XX('/Users/ahoffmann/gyacomo/results/paper_3/DTT_rho85/3x2x192x48x32');
half = ceil(numel(t_PT)/2);
Gi_avg_PT = mean(Gi_PT(half:end));
Ge_avg_PT = mean(Ge_PT(half:end));
Qi_avg_PT = mean(Qi_PT(half:end));
Qe_avg_PT = mean(Qe_PT(half:end));


[t_NT, Gi_NT, Qi_NT, Ge_NT, Qe_NT] = read_flux_out_XX('/Users/ahoffmann/gyacomo/results/paper_3/DTT_rho85/3x2x192x48x32_NT');
half = ceil(numel(t_NT)/2);
Gi_avg_NT = mean(Gi_NT(half:end));
Ge_avg_NT = mean(Ge_NT(half:end));
Qi_avg_NT = mean(Qi_NT(half:end));
Qe_avg_NT = mean(Qe_NT(half:end));

prct_change = @(x,y) abs(y-x)/abs(x)*100;
nmvm = 10;
plt_ = @(t,f,stl,nme) plot(movmean(t,nmvm),movmean(f,nmvm), stl,'DisplayName',nme);

figure

nr = 1;
if 1
    nr = 2;
    subplot(223)
    plt_(t_PT,Gi_PT, '-r','ions PT'); hold on;
    plt_(t_NT,Gi_NT,'--r','ions NT'); hold on;
    xlabel('$tc_s/R$'); ylabel('$\Gamma_{x}$');
    legend('show')
    title(['Ion particle flux, NT/NP=',sprintf('%3.1f',prct_change(Gi_avg_PT,Gi_avg_NT)),'$\%$'])
    subplot(224)
    plt_(t_PT,Ge_PT, '-b','electrons PT'); hold on;
    plt_(t_NT,Ge_NT,'--b','electrons NT'); hold on;
    xlabel('$tc_s/R$'); ylabel('$\Gamma_{x}$');
    legend('show')
    title(['Elec. particle flux, NT/NP=',sprintf('%3.1f',prct_change(Ge_avg_PT,Ge_avg_NT)),'$\%$'])
end

subplot(nr,2,1)
plt_(t_PT,Qi_PT, '-r','ions PT'); hold on;
plt_(t_NT,Qi_NT,'--r','ions NT'); hold on;
xlabel('$tc_s/R$'); ylabel('$Q_{x}$');
legend('show')
title(['Ion heat flux, NT/NP=',sprintf('%3.1f',prct_change(Qi_avg_PT,Qi_avg_NT)),'$\%$'])

subplot(nr,2,2)
plt_(t_PT,Qe_PT, '-b','electrons PT'); hold on;
plt_(t_NT,Qe_NT,'--b','electrons NT'); hold on;
xlabel('$tc_s/R$'); ylabel('$Q_{x}$');
legend('show')
title(['Elec. heat flux, NT/NP=',sprintf('%3.1f',prct_change(Qe_avg_PT,Qe_avg_NT)),'$\%$'])

disp(['PT total heat flux : ',num2str(Qi_avg_PT + Qe_avg_PT)]);
disp(['NT total heat flux : ',num2str(Qi_avg_NT + Qe_avg_NT)]);
