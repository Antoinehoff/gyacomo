function [ FIGURE ] = plot_fa( DATA, OPTIONS )

FIGURE.fig = figure; FIGURE.FIGNAME = ['f_a_',DATA.PARAMS];
switch OPTIONS.Z
    case 'avg'
        zcomp = ' z-avg';
    otherwise
        zcomp = [' z=',sprintf('%2.2f',DATA.z(OPTIONS.Z))];
end
if OPTIONS.ONED
    [s,x,fsa,fxa] = compute_fa_1D(DATA, OPTIONS); 
    [~,it] = min(abs(OPTIONS.T-DATA.Ts5D)); 
    subplot(1,2,1)
        plot(s,fsa); hold on
        legend(OPTIONS.SPECIE)
        xlabel('$v_\parallel, (\mu=0)$'); ylabel('$\langle |f_a|^2\rangle_{xy}^{1/2}$'); 
        title(DATA.param_title); 
    subplot(1,2,2)
        plot(x,fxa); hold on;
        legend(OPTIONS.SPECIE)
        xlabel('$\mu, (v_\parallel=0)$'); ylabel('$\langle |f_a|^2\rangle_{xy}^{1/2}$'); 
        if numel(it) == 1
            title(['t=',num2str(DATA.Ts5D(it))]);
        else
            title(['average $t\in$[',num2str(DATA.Ts5D(min(it))),',',num2str(DATA.Ts5D(max(it))),']']);
        end

else
    [SS,XX,FFa] = compute_fa_2D(DATA, OPTIONS); 
    [~,it] = min(abs(OPTIONS.T-DATA.Ts5D)); 
        switch OPTIONS.PLT_FCT
            case 'contour'
            contour(SS,XX,FFa',128);
            case 'pcolor'
            pclr = pcolor(SS,XX,FFa'); set(pclr, 'edgecolor','none'); shading interp
        end
        xlabel('$v_\parallel$'); ylabel('$\mu$');
        legend(['$\langle |f_',OPTIONS.SPECIE,'|^2\rangle_{xy}^{1/2}$',zcomp])
        if numel(it) == 1
            title(['HeLaZ''$\langle |f_',OPTIONS.SPECIE...
                ,'|^2\rangle_{xy}^{1/2}$',zcomp...
                ,', $t=$[',sprintf('%1.1f',DATA.Ts5D(it)),']']);
        else
            title(['HeLaZ''$\langle |f_',OPTIONS.SPECIE,...
                '|^2\rangle_{xy}^{1/2}$',zcomp,...
                ', average $t\in$[',sprintf('%1.1f',DATA.Ts5D(min(it))),',',sprintf('%1.1f',DATA.Ts5D(max(it))),']']);
        end
end

end

