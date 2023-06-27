function [FIG] = plot_ballooning(data,options)
    FIG.fig = figure; iplot = 1;
    [~,it0] = min(abs(data.Ts3D - options.time_2_plot(1)));
    [~,it1] = min(abs(data.Ts3D - options.time_2_plot(end)));
    [~,ikyarray] = min(abs(data.grids.ky - options.kymodes));
    ikyarray = unique(ikyarray);
    dz = data.grids.z(2) - data.grids.z(1);
    phi_real=real(data.PHI(:,:,:,it1));
    phi_imag=imag(data.PHI(:,:,:,it1));
    ncol = 1;
    if data.inputs.BETA > 0
        psi_real=real(data.PSI(:,:,:,it1));
        psi_imag=imag(data.PSI(:,:,:,it1)); 
        ncol = 2;
    end
    if options.PLOT_KP
        ncol = 3;
    end
    % % Apply ballooning transform
    % if(data.grids.Nkx > 1)
    %     nexc = round(data.grids.ky(2)*data.inputs.SHEAR*2*pi/data.grids.kx(2));
    % else
    %     nexc = 1;
    % end
    for iky=ikyarray

        [phib_real,b_angle] = ballooning_representation(phi_real,data.inputs.SHEAR,data.grids.Npol,data.grids.kx,iky,data.grids.ky,data.grids.z);
        [phib_imag,~]       = ballooning_representation(phi_imag,data.inputs.SHEAR,data.grids.Npol,data.grids.kx,iky,data.grids.ky,data.grids.z);
        phib = phib_real(:) + 1i * phib_imag(:);
        % normalize real and imaginary parts at chi =0
        if options.normalized
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
        plot(b_angle/pi, phib_real_norm,'-b'); hold on;
        plot(b_angle/pi, phib_imag_norm ,'-r');
        plot(b_angle/pi, sqrt(phib_real_norm .^2 + phib_imag_norm.^2),'-k');
        legend('real','imag','norm')
        xlabel('$\chi / \pi$')
        ylabel('$\phi_B (\chi)$');
        title(['$k_y=',sprintf('%2.2f',data.grids.ky(iky)),...
               ',t_{avg}\in [',sprintf('%1.1f',data.Ts3D(it0)),',',sprintf('%1.1f',data.Ts3D(it1)),']$']);
        % z domain start end point
        for i_ = 2*[1:data.grids.Nkx-1]-(data.grids.Nkx+1)
            xline(data.grids.Npol*(i_),'HandleVisibility','off'); hold on;
        end
        for i_ = 2*[2:data.grids.Nkx]-(data.grids.Nkx+1)
            xline(data.grids.Npol*(i_)-dz/pi,'HandleVisibility','off'); hold on;
        end
        if data.inputs.BETA > 0         
            [psib_real,b_angle] = ballooning_representation(psi_real,data.inputs.SHEAR,data.grids.Npol,data.grids.kx,iky,data.grids.ky,data.grids.z);
            [psib_imag,~]       = ballooning_representation(psi_imag,data.inputs.SHEAR,data.grids.Npol,data.grids.kx,iky,data.grids.ky,data.grids.z);

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
            for i_ = 2*[1:data.grids.Nkx-1]-(data.grids.Nkx+1)
                xline(data.grids.Npol*(i_),'HandleVisibility','off'); hold on;
            end
            for i_ = 2*[2:data.grids.Nkx]-(data.grids.Nkx+1)
                xline(data.grids.Npol*(i_)-dz/pi,'HandleVisibility','off'); hold on;
            end
            xlabel('$\chi / \pi$')
            ylabel('$\psi_B (\chi)$');
            title(['$k_y=',sprintf('%2.2f',data.grids.ky(iky)),...
                   ',t_{avg}\in [',sprintf('%1.1f',data.Ts3D(it0)),',',sprintf('%1.1f',data.Ts3D(it1)),']$']);
        end
        
        if options.PLOT_KP
            kperp      = h5read(data.outfilenames{1},'/data/grid/coordkp');
            [kperpb,b_angle] = ballooning_representation(kperp,data.inputs.SHEAR,data.grids.Npol,data.grids.kx,iky,data.grids.ky,data.grids.z);
            subplot(numel(ikyarray),ncol,ncol*(iplot-1)+3)
            plot(b_angle/pi, kperpb,'-k'); hold on;
            % z domain start end point
            for i_ = 2*[1:data.grids.Nkx-1]-(data.grids.Nkx+1)
                xline(data.grids.Npol*(i_),'HandleVisibility','off'); hold on;
            end
            for i_ = 2*[2:data.grids.Nkx]-(data.grids.Nkx+1)
                xline(data.grids.Npol*(i_)-dz/pi,'HandleVisibility','off'); hold on;
            end
            xlabel('$\chi / \pi$')
            ylabel('$k_\perp (\chi)$');
            title(['$k_y=',sprintf('%2.2f',data.grids.ky(iky)),...
                   ',t_{avg}\in [',sprintf('%1.1f',data.Ts3D(it0)),',',sprintf('%1.1f',data.Ts3D(it1)),']$']);

        end
        iplot = iplot + 1;
    end
end
