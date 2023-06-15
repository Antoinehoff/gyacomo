function [FIG] = plot_ballooning(data,options)
    FIG.fig = figure; iplot = 1;
    [~,it0] = min(abs(data.Ts3D - options.time_2_plot(1)));
    [~,it1] = min(abs(data.Ts3D - options.time_2_plot(end)));
    [~,ikyarray] = min(abs(data.grids.ky - options.kymodes));
    ikyarray = unique(ikyarray);
    phi_real=real(data.PHI(:,:,:,it1));
    phi_imag=imag(data.PHI(:,:,:,it1));
    ncol = 1;
    if data.inputs.BETA > 0
        psi_real=real(data.PSI(:,:,:,it1));
        psi_imag=imag(data.PSI(:,:,:,it1)); 
        ncol = 2;
    end
    % Apply ballooning transform
    if(data.grids.Nkx > 1)
        nexc = round(data.grids.ky(2)*data.inputs.SHEAR*2*pi/data.grids.kx(2));
    else
        nexc = 1;
    end
    for iky=ikyarray
        dims = size(phi_real);
        Nkx  = dims(2);
        is   = max(1,iky-1);
        Npi  = (Nkx-1)-2*nexc*(is-1);
        if(Npi <= 1)
            ordered_ikx = 1;
        elseif(mod(Nkx,2) == 0)
            tmp_ = (Nkx-is+1):-is:(Nkx/2+2);
            ordered_ikx = [tmp_(end:-1:1), 1:is:(Nkx/2)];
        else
            Np_ = (Nkx+1)/(2*is);
            ordered_ikx = [(Np_+1):Nkx 1:Np_];
        end
        try
            idx=data.grids.kx./data.grids.kx(2);
        catch
            idx=0;
        end
        idx= idx(ordered_ikx);
        Nkx = numel(idx);

        phib_real = zeros(  Nkx*dims(3)  ,1);
        phib_imag = phib_real;
        b_angle   = phib_real;

        for i_ =1:Nkx
            start_ =  (i_-1)*dims(3) +1;
            end_ =  i_*dims(3);
            ikx  = ordered_ikx(i_);
            phib_real(start_:end_)  = phi_real(iky,ikx,:);
            phib_imag(start_:end_)  = phi_imag(iky,ikx,:);
        end

        % Define ballooning angle
        coordz = data.grids.z;
        for i_ =1: Nkx
            for iz=1:dims(3)
                ii = dims(3)*(i_-1) + iz;
                b_angle(ii) =coordz(iz) + 2*pi*data.grids.Npol*idx(i_)./is;
            end
        end
        
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
        for i_ = 2*[1:data.grids.Nkx]-(data.grids.Nkx+1)
          xline(data.grids.Npol*(i_),'HandleVisibility','off'); hold on;
        end
        if data.inputs.BETA > 0         
            psib_real = zeros(  Nkx*dims(3)  ,1);
            psib_imag = psib_real;
            for i_ =1:Nkx
                start_ =  (i_-1)*dims(3) +1;
                end_ =  i_*dims(3);
                ikx  = ordered_ikx(i_);
                psib_real(start_:end_)  = psi_real(iky,ikx,:);
                psib_imag(start_:end_)  = psi_imag(iky,ikx,:);
            end
            psib = psib_real(:) + 1i * psib_imag(:);
            psib_norm = psib / normalization;
            psib_real_norm  = real( psib_norm);
            psib_imag_norm  = imag( psib_norm);   
            subplot(numel(ikyarray),ncol,ncol*(iplot-1)+2)
            plot(b_angle/pi, psib_real_norm,'-b'); hold on;
            plot(b_angle/pi, psib_imag_norm ,'-r');
            plot(b_angle/pi, abs(psib_norm),'-k');
            legend('real','imag','norm')
            for i_ = 2*[1:data.grids.Nkx]-(data.grids.Nkx+1)
                xline(data.grids.Npol*(i_),'HandleVisibility','off'); hold on;
            end
            xlabel('$\chi / \pi$')
            ylabel('$\psi_B (\chi)$');
            title(['$k_y=',sprintf('%2.2f',data.grids.ky(iky)),...
                   ',t_{avg}\in [',sprintf('%1.1f',data.Ts3D(it0)),',',sprintf('%1.1f',data.Ts3D(it1)),']$']);
        end
        
        iplot = iplot + 1;
    end
end
