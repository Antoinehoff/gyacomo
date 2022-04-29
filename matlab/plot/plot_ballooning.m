function [FIG] = plot_ballooning(data,options)
    
    [~,it] = min(abs(data.Ts3D - options.time_2_plot));
    phi_real=(real(data.PHI(:,:,:,it)));
    phi_imag=(imag(data.PHI(:,:,:,it)));
    % Apply baollooning tranform
    for iky=2
        dims = size(phi_real);
        phib_real = zeros(  dims(1)*dims(3)  ,1);
        phib_imag= phib_real;
        b_angle = phib_real;

        midpoint = floor((dims(1)*dims(3) )/2)+1;
        for ip =1: dims(1)
            start_ =  (ip -1)*dims(3) +1;
            end_ =  ip*dims(3);
            phib_real(start_:end_)  = phi_real(ip,iky,:);
            phib_imag(start_:end_)  = phi_imag(ip,iky, :);
        end

        % Define ballooning angle
        Nkx = numel(data.kx)-1; coordz = data.z;
        idx = -Nkx:1:Nkx;
        for ip =1: dims(1)
            for iz=1:dims(3)
                ii = dims(3)*(ip -1) + iz;
                b_angle(ii) =coordz(iz) + 2*pi*idx(ip);
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


        FIG.fig = figure; hold all;
        plot(b_angle/pi, phib_real_norm,'-b');
        plot(b_angle/pi, phib_imag_norm ,'-r');
        plot(b_angle/pi, sqrt(phib_real_norm .^2 + phib_imag_norm.^2),'-k');
        legend('real','imag','norm')
        xlabel('$\chi / \pi$')
        ylabel('$\phi_B (\chi)$');
        title(['HeLaZ ballooning, t=',num2str(data.Ts3D(it))]);
    end
end
