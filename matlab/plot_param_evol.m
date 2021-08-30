figure; hold on; set(gcf, 'Position',  [100, 100, 400, 1500])
subplot(4,1,1);
    plot(TJOB_SE,NU_EVOL,'DisplayName','$\nu$');
    xlim([TJOB_SE(1) TJOB_SE(end)]);legend('show');grid on;
    
subplot(4,1,2);
    plot(TJOB_SE,MU_EVOL,'DisplayName','$\mu$');
    xlim([TJOB_SE(1) TJOB_SE(end)]);legend('show');grid on;
    
subplot(4,1,3);
    plot(TJOB_SE,1./ETAN_EVOL,'DisplayName','$\eta$');
    xlim([TJOB_SE(1) TJOB_SE(end)]);legend('show');grid on;

subplot(4,1,4);
    plot(TJOB_SE,L_EVOL,'DisplayName','$L$');
    xlim([TJOB_SE(1) TJOB_SE(end)]);legend('show');grid on;

 xlabel('$t c_s/R$');
