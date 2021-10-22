fig = figure; FIGNAME = 'linear_study';
hold on; set(gcf, 'Position',  [100, 100, 400, 1500])
subplot(4,1,1);
    title('Parameter evolution'); hold on;
    yyaxis left
    plot(TJOB_SE,NU_EVOL,'DisplayName','$\nu$');
    yyaxis right
    plot(TJOB_SE,CO_EVOL,'DisplayName','CO');
    xlim([TJOB_SE(1) TJOB_SE(end)]);legend('show');grid on;
    xticks([]);
    plot_tjob_lines(TJOB_SE,ylim)
subplot(4,1,2);
    plot(TJOB_SE,MU_EVOL,'DisplayName','$\mu$'); hold on;
    xlim([TJOB_SE(1) TJOB_SE(end)]);legend('show');grid on;
    xticks([]);
    plot_tjob_lines(TJOB_SE,ylim)

subplot(4,1,3);
    yyaxis left
    plot(TJOB_SE,K_N_EVOL,'DisplayName','$\nabla n$'); hold on;
    yyaxis right
    plot(TJOB_SE,K_T_EVOL,'DisplayName','$\nabla T$'); hold on;
    xlim([TJOB_SE(1) TJOB_SE(end)]);legend('show');grid on;
    xticks([]);
    plot_tjob_lines(TJOB_SE,ylim)

subplot(4,1,4);
    plot(TJOB_SE,L_EVOL,'DisplayName','$L$'); hold on;
    xlim([TJOB_SE(1) TJOB_SE(end)]);legend('show');grid on;
    plot_tjob_lines(TJOB_SE,ylim)

 xlabel('$t c_s/R$');
 saveas(fig,[BASIC.RESDIR,'param_evol.png']);

 function [] = plot_tjob_lines(TJOB,limits)
     for i = 2:numel(TJOB)-1
        plot(TJOB(i)*[1 1],limits,'--k')
     end
 end
