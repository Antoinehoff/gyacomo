function [fig] = plot_param_evol(data)
TJOB_SE  = data.TJOB_SE;
NU_EVOL  = data.NU_EVOL;
DT_EVOL  = data.DT_EVOL;
CO_EVOL  = data.CO_EVOL;
MUx_EVOL  = data.MUx_EVOL;
MUy_EVOL  = data.MUy_EVOL;
MUz_EVOL  = data.MUz_EVOL;
% K_T_EVOL = data.K_T_EVOL;
K_N_EVOL = data.K_N_EVOL;
L_EVOL   = data.L_EVOL;

CO_LIST = {};
for i_ = 1:numel(CO_EVOL)/6
    CO_LIST{i_} = CO_EVOL(6*(i_-1)+1:6*i_);
end
fig = figure;
hold on; set(gcf, 'Position',  [100, 100, 400, 1500])
subplot(4,1,1);
    title('Parameter evolution'); hold on;
    yyaxis left
    plot(TJOB_SE,NU_EVOL,'DisplayName','$\nu$');
    for i_ = 1:numel(CO_EVOL)/12
        k_ = 2*(i_-1)+1;
    text(double(TJOB_SE(k_)),NU_EVOL(k_)*2.0,CO_LIST{k_});
    end
    yyaxis right
    plot(TJOB_SE,DT_EVOL,'DisplayName','dt');
    xlim([TJOB_SE(1) TJOB_SE(end)]);legend('show');grid on;
    xticks([]);
    plot_tjob_lines(TJOB_SE,ylim)
subplot(4,1,2);
    plot(TJOB_SE,MUx_EVOL,'DisplayName','$\mu_x$'); hold on;
    plot(TJOB_SE,MUy_EVOL,'DisplayName','$\mu_y$'); hold on;
    plot(TJOB_SE,MUz_EVOL,'DisplayName','$\mu_z$'); hold on;
    xlim([TJOB_SE(1) TJOB_SE(end)]);legend('show');grid on;
    xticks([]);
    plot_tjob_lines(TJOB_SE,ylim)

subplot(4,1,3);
    yyaxis left
    plot(TJOB_SE,K_N_EVOL,'DisplayName','$\nabla n$'); hold on;
    yyaxis right
%     plot(TJOB_SE,K_T_EVOL,'DisplayName','$\nabla T$'); hold on;
    xlim([TJOB_SE(1) TJOB_SE(end)]);legend('show');grid on;
    xticks([]);
    plot_tjob_lines(TJOB_SE,ylim)

subplot(4,1,4);
    plot(TJOB_SE,L_EVOL,'DisplayName','$L$'); hold on;
    xlim([TJOB_SE(1) TJOB_SE(end)]);legend('show');grid on;
    plot_tjob_lines(TJOB_SE,ylim)

 xlabel('$t c_s/R$');

 function [] = plot_tjob_lines(TJOB,limits)
     for i = 2:numel(TJOB)-1
        plot(TJOB(i)*[1 1],limits,'--k')
     end
 end

end