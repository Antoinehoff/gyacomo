function [ FIGURE ] = plot_fa( DATA, OPTIONS )

FIGURE.fig = figure; FIGURE.FIGNAME = ['f_a_',DATA.PARAMS];

if OPTIONS.ONED
    OPTIONS.SPECIE = 'i';
    [~,~,fsi,fxi] = compute_fa_1D(DATA, OPTIONS); 
    OPTIONS.SPECIE = 'e';
    [s,x,fse,fxe] = compute_fa_1D(DATA, OPTIONS); 
    [~,it] = min(abs(OPTIONS.T-DATA.Ts5D)); 
    subplot(1,2,1)
        plot(s,fse); hold on
        plot(s,fsi);
        legend('e','i')
        xlabel('$v_\parallel, (\mu=0)$'); ylabel('$\langle |f_a|^2\rangle_{xy}^{1/2}$'); 
        title(DATA.param_title); 
    subplot(1,2,2)
        plot(x,fxe); hold on;
        plot(x,fxi);
        legend('e','i')
        xlabel('$\mu, (v_\parallel=0)$'); ylabel('$\langle |f_a|^2\rangle_{xy}^{1/2}$'); 
        if numel(it) == 1
            title(['t=',num2str(DATA.Ts5D(it))]);
        else
            title(['average $t\in$[',num2str(DATA.Ts5D(min(it))),',',num2str(DATA.Ts5D(max(it))),']']);
        end

else
    OPTIONS.SPECIE = 'i';
    [~,~,FFi] = compute_fa_2D(DATA, OPTIONS); 
    OPTIONS.SPECIE = 'e';
    [SS,XX,FFe] = compute_fa_2D(DATA, OPTIONS); 
    [~,it] = min(abs(OPTIONS.T-DATA.Ts5D)); 
    subplot(1,2,1)
        if OPTIONS.CTR
        contour(SS,XX,FFi',128);
        else
        pclr = pcolor(SS,XX,FFi'); set(pclr, 'edgecolor','none'); shading interp
        end
        xlabel('$v_\parallel$'); ylabel('$\mu$');
        legend('$\langle |f_i|^2\rangle_{xy}^{1/2}$')
        title(DATA.param_title); 
        % FF = log10(FF);
    subplot(1,2,2)
        if OPTIONS.CTR
        contour(SS,XX,FFe',128);
        else
        pclr = pcolor(SS,XX,FFe'); set(pclr, 'edgecolor','none'); shading interp
        end
        legend('$\langle |f_e|^2\rangle_{xy}^{1/2}$')
        xlabel('$v_\parallel$'); ylabel('$\mu$');
        if numel(it) == 1
            title(['t=',num2str(DATA.Ts5D(it))]);
        else
            title(['average $t\in$[',num2str(DATA.Ts5D(min(it))),',',num2str(DATA.Ts5D(max(it))),']']);
        end
end

end

