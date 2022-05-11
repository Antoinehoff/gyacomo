function [FIG] = plot_ballooning(data,options)
    FIG.fig = figure; iplot = 1;
    [~,it] = min(abs(data.Ts3D - options.time_2_plot));
    [~,ikyarray] = min(abs(data.ky - options.kymodes));
    phi_real=(real(data.PHI(:,:,:,it)));
    phi_imag=(imag(data.PHI(:,:,:,it)));
    % Apply baollooning tranform
    for iky=ikyarray
        dims = size(phi_real);
        phib_real = zeros(  dims(2)*dims(3)  ,1);
        phib_imag= phib_real;
        b_angle = phib_real;

        for ikx =1: dims(2)
            start_ =  (ikx -1)*dims(3) +1;
            end_ =  ikx*dims(3);
            phib_real(start_:end_)  = phi_real(iky,ikx,:);
            phib_imag(start_:end_)  = phi_imag(iky,ikx,:);
        end

        % Define ballooning angle
        Nkx = numel(data.kx)-1; coordz = data.z;
        idx = -Nkx:1:Nkx;
        for ikx =1: dims(2)
            for iz=1:dims(3)
                ii = dims(3)*(ikx -1) + iz;
                b_angle(ii) =coordz(iz) + 2*pi*idx(ikx);
            end
        end
        
        phib = phib_real(:) + 1i * phib_imag(:);
        % normalize real and imaginary parts at chi =0
        if options.normalized
            [~,idxLFS] = min(abs(b_angle -0));
            normalization = phib( idxLFS);
        else
            normalization = 1;
        end
        phib_norm = phib / normalization    ;
        phib_real_norm  = real( phib_norm);%phib_real(:)/phib_real(idxLFS);
        phib_imag_norm  = imag( phib_norm);%phib_imag(:)/ phib_imag(idxLFS);

        subplot(numel(ikyarray),1,iplot)
        plot(b_angle/pi, phib_real_norm,'-b'); hold on;
        plot(b_angle/pi, phib_imag_norm ,'-r');
        plot(b_angle/pi, sqrt(phib_real_norm .^2 + phib_imag_norm.^2),'-k');
        legend('real','imag','norm')
        xlabel('$\chi / \pi$')
        ylabel('$\phi_B (\chi)$');
        title(['ky=',sprintf('%1.1f',data.ky(iky)),...
               ',t=',sprintf('%1.1f',data.Ts3D(it))]);
        iplot = iplot + 1;
    end
end
