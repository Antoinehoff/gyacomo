function [ FIGURE ] = plot_fa( DATA, OPTIONS )

FIGURE.fig = figure; FIGURE.FIGNAME = ['f_a_',DATA.PARAMS];
switch OPTIONS.iz
    case 'avg'
        zcomp = ' z-avg';
    otherwise
        zcomp = [' z=',sprintf('%2.2f',DATA.z(OPTIONS.iz))];
end
if OPTIONS.ONED
    [s,x,fsa,fxa] = compute_fa_1D(DATA, OPTIONS); 
    [~,it] = min(abs(OPTIONS.T-DATA.Ts5D)); 
    FIGURE.ax1 = subplot(1,2,1,'parent',FIGURE.fig);
        plot(s,fsa); hold on
        legend(OPTIONS.SPECIES)
        xlabel('$v_\parallel, (\mu=0)$'); ylabel(['$\langle |f_a|^2\rangle_{xy}^{1/2}$, ',zcomp]); 
        title(DATA.param_title); 
    FIGURE.ax2 = subplot(1,2,2,'parent',FIGURE.fig);
        plot(x,fxa); hold on;
        legend(OPTIONS.SPECIES)
        xlabel('$\mu, (v_\parallel=0)$'); ylabel(['$\langle |f_a|^2\rangle_{xy}^{1/2}$, ',zcomp]);
        if numel(it) == 1
            title(['t=',num2str(DATA.Ts5D(it))]);
        else
            title(['average $t\in$[',num2str(DATA.Ts5D(min(it))),',',num2str(DATA.Ts5D(max(it))),']']);
        end

else
    FIGURE.ax1 = subplot(1,1,1,'parent',FIGURE.fig);
    [SS,XX,FFa] = compute_fa_2D(DATA, OPTIONS);  sz = size(SS); FFa = FFa';
    [~,it] = min(abs(OPTIONS.T-DATA.Ts5D)); 
        switch OPTIONS.PLT_FCT
            case 'contour'
            contour(SS,XX,FFa,sum(sz)/2);
            xlabel('$v_\parallel$'); ylabel('$\mu$');
            case 'pcolor'
            pclr = pcolor(SS,XX,FFa); set(pclr, 'edgecolor','none'); shading interp
            xlabel('$v_\parallel$'); ylabel('$\mu$');
            case 'contourf'
            contourf(SS,XX,FFa,sum(sz)/2);    
            xlabel('$v_\parallel$'); ylabel('$\mu$');
            case 'surf'
            surf(SS,XX,FFa); 
            xlabel('$v_\parallel$'); ylabel('$\mu$');
            case 'surfvv'
            surf([SS(end:-1:1,:) SS ],[-sqrt(XX(end:-1:1,:)) sqrt(XX)],[FFa(end:-1:1,:) FFa]);
            xlabel('$v_\parallel$'); ylabel('$v_\perp$');
            xlim([min(OPTIONS.SPAR) max(OPTIONS.SPAR)]);
            ylim(sqrt(max(OPTIONS.XPERP))*[-1 1]);
        end
        legend(['$\langle |f_',OPTIONS.SPECIES,'|^2\rangle_{xy}^{1/2}$',zcomp])
        if numel(it) == 1
            title(['HeLaZ''$\langle |f_',OPTIONS.SPECIES...
                ,'|^2\rangle_{xy}^{1/2}$',zcomp...
                ,', $t=$[',sprintf('%1.1f',DATA.Ts5D(it)),']']);
        else
            title(['HeLaZ''$\langle |f_',OPTIONS.SPECIES,...
                '|^2\rangle_{xy}^{1/2}$',zcomp,...
                ', average $t\in$[',sprintf('%1.1f',DATA.Ts5D(min(it))),',',sprintf('%1.1f',DATA.Ts5D(max(it))),']']);
        end
end

end

