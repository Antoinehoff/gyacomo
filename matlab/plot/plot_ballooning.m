function [FIG, b_angle, phib, psib, ky_] = plot_ballooning(DATA,OPTIONS)
    FIG.fig = figure; iplot = 1;
    [~,it0] = min(abs(DATA.Ts3D - OPTIONS.time_2_plot(1)));
    [~,it1] = min(abs(DATA.Ts3D - OPTIONS.time_2_plot(end)));
    [~,ikyarray] = min(abs(DATA.grids.ky - OPTIONS.kymodes));
    ikyarray = unique(ikyarray);
    dz = DATA.grids.z(2) - DATA.grids.z(1);
    phi_real=real(DATA.PHI(:,:,:,it1));
    phi_imag=imag(DATA.PHI(:,:,:,it1));
    ncol = 1;
    if DATA.inputs.BETA > 0
        psi_real=real(DATA.PSI(:,:,:,it1));
        psi_imag=imag(DATA.PSI(:,:,:,it1)); 
        ncol = 2;
    end
    if OPTIONS.PLOT_KP
        ncol = 3;
    end
    % % Apply ballooning transform
    % if(DATA.grids.Nkx > 1)
    %     nexc = round(DATA.grids.ky(2)*DATA.inputs.SHEAR*2*pi/DATA.grids.kx(2));
    % else
    %     nexc = 1;
    % end
    nline = 1;
    for iky=ikyarray

        [phib_real,b_angle] = ballooning_representation(phi_real,DATA.inputs.SHEAR,DATA.grids.Npol,DATA.grids.kx,iky,DATA.grids.ky,DATA.grids.z);
        [phib_imag,~]       = ballooning_representation(phi_imag,DATA.inputs.SHEAR,DATA.grids.Npol,DATA.grids.kx,iky,DATA.grids.ky,DATA.grids.z);
        phib = phib_real(:) + 1i * phib_imag(:);
        % normalize real and imaginary parts at chi =0
        if OPTIONS.normalized
            [~,idxLFS] = min(abs(b_angle -0));
            normalization = (phib( idxLFS));
        else
            normalization = 1;
        end
        phib_norm = phib / normalization;
        phib_real_norm  = real( phib_norm);
        phib_imag_norm  = imag( phib_norm);
        %
        subplot(numel(ikyarray),ncol,ncol*(iplot-1)+1)
        plot(b_angle/pi, phib_real_norm,'-b','DisplayName','$|\phi_B (\chi)|$'); hold on;
        plot(b_angle/pi, phib_imag_norm ,'-r');
        plot(b_angle/pi, sqrt(phib_real_norm .^2 + phib_imag_norm.^2),'-k');
        legend('real','imag','norm')
        xlabel('$\chi / \pi$')
        ylabel('$\phi_B (\chi)$');
        title(['$k_y=',sprintf('%2.2f',DATA.grids.ky(iky)),...
               ',t_{avg}\in [',sprintf('%1.1f',DATA.Ts3D(it0)),',',sprintf('%1.1f',DATA.Ts3D(it1)),']$']);
        % z domain start end point
        for i_ = 2*[1:DATA.grids.Nkx-1]-(DATA.grids.Nkx+1)
            xline(DATA.grids.Npol*(i_),'HandleVisibility','off'); hold on;
        end
        for i_ = 2*[2:DATA.grids.Nkx]-(DATA.grids.Nkx+1)
            xline(DATA.grids.Npol*(i_)-dz/pi,'HandleVisibility','off'); hold on;
        end
        if DATA.inputs.BETA > 0         
            [psib_real,b_angle] = ballooning_representation(psi_real,DATA.inputs.SHEAR,DATA.grids.Npol,DATA.grids.kx,iky,DATA.grids.ky,DATA.grids.z);
            [psib_imag,~]       = ballooning_representation(psi_imag,DATA.inputs.SHEAR,DATA.grids.Npol,DATA.grids.kx,iky,DATA.grids.ky,DATA.grids.z);

            psib = psib_real(:) + 1i * psib_imag(:);
            psib_norm = psib / normalization;
            psib_real_norm  = real( psib_norm);
            psib_imag_norm  = imag( psib_norm);   
            subplot(numel(ikyarray),ncol,ncol*(iplot-1)+2)
            plot(b_angle/pi, psib_real_norm,'-b'); hold on;
            plot(b_angle/pi, psib_imag_norm ,'-r');
            plot(b_angle/pi, abs(psib_norm),'-k');

            legend('real','imag','norm')
            % z domain start end point
            for i_ = 2*[1:DATA.grids.Nkx-1]-(DATA.grids.Nkx+1)
                xline(DATA.grids.Npol*(i_),'HandleVisibility','off'); hold on;
            end
            for i_ = 2*[2:DATA.grids.Nkx]-(DATA.grids.Nkx+1)
                xline(DATA.grids.Npol*(i_)-dz/pi,'HandleVisibility','off'); hold on;
            end
            xlabel('$\chi / \pi$')
            ylabel('$\psi_B (\chi)$');
            title(['$k_y=',sprintf('%2.2f',DATA.grids.ky(iky)),...
                   ',t_{avg}\in [',sprintf('%1.1f',DATA.Ts3D(it0)),',',sprintf('%1.1f',DATA.Ts3D(it1)),']$']);
        end
        
        if OPTIONS.PLOT_KP
            kperp      = h5read(DATA.outfilenames{1},'/DATA/grid/coordkp');
            [kperpb,b_angle] = ballooning_representation(kperp,DATA.inputs.SHEAR,DATA.grids.Npol,DATA.grids.kx,iky,DATA.grids.ky,DATA.grids.z);
            subplot(numel(ikyarray),ncol,ncol*(iplot-1)+3)
            plot(b_angle/pi, kperpb,'-k'); hold on;
            % z domain start end point
            for i_ = 2*[1:DATA.grids.Nkx-1]-(DATA.grids.Nkx+1)
                xline(DATA.grids.Npol*(i_),'HandleVisibility','off'); hold on;
            end
            for i_ = 2*[2:DATA.grids.Nkx]-(DATA.grids.Nkx+1)
                xline(DATA.grids.Npol*(i_)-dz/pi,'HandleVisibility','off'); hold on;
            end
            xlabel('$\chi / \pi$')
            ylabel('$k_\perp (\chi)$');
            title(['$k_y=',sprintf('%2.2f',DATA.grids.ky(iky)),...
                   ',t_{avg}\in [',sprintf('%1.1f',DATA.Ts3D(it0)),',',sprintf('%1.1f',DATA.Ts3D(it1)),']$']);

        end
        iplot = iplot + 1;
    end
    if ~OPTIONS.SHOWFIG
        close
    end
    ky_ = DATA.grids.ky(iky);
end

function [plots, ncurve] = add_to_plots(plots,plotname,x,xname,y,yname,ncurve)
            plots.(plotname).x = x;
            plots.(plotname).y = y;
            plots.(plotname).xname = xname;
            plots.(plotname).yname = yname;
end