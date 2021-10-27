[~,it] = min(abs(Ts3D - time_2_plot));
phi_real=(real(PHI(:,:,:,it)));
phi_imag=(imag(PHI(:,:,:,it)));
% Apply baollooning tranform
for iky=2
    dims = size(phi_real);
    phib_real = zeros(  dims(1)*Nz  ,1);
    phib_imag= phib_real;
    b_angle = phib_real;
    
    midpoint = floor((dims(1)*Nz )/2)+1;
    for ip =1: dims(1)
        start_ =  (ip -1)*Nz +1;
        end_ =  ip*Nz;
        phib_real(start_:end_)  = phi_real(ip,iky,:);
        phib_imag(start_:end_)  = phi_imag(ip,iky, :);
    end
    
    % Define ballooning angle
    Nkx = numel(kx)-1; coordz = z;
    idx = -Nkx:1:Nkx;
    for ip =1: dims(1)
        for iz=1:Nz
            ii = Nz*(ip -1) + iz;
            b_angle(ii) =coordz(iz) + 2*pi*idx(ip);
        end
    end
    
    % normalize real and imaginary parts at chi =0
    [~,idxLFS] = min(abs(b_angle -0));
    phib = phib_real(:) + 1i * phib_imag(:);
    % Normalize to the outermid plane
    phib_norm = phib(:);% / phib( idxLFS)    ;
    phib_real_norm(:)  = real( phib_norm);%phib_real(:)/phib_real(idxLFS);
    phib_imag_norm(:)  = imag( phib_norm);%phib_imag(:)/ phib_imag(idxLFS);
    
    
    figure; hold all;
    plot(b_angle/pi, phib_real_norm,'-b');
    plot(b_angle/pi, phib_imag_norm ,'-r');
    plot(b_angle/pi, sqrt(phib_real_norm .^2 + phib_imag_norm.^2),'-k');
    legend('real','imag','norm')
    xlabel('$\chi / \pi$')
    ylabel('$\phi_B (\chi)$');
    title(['HeLaZ(-) molix(o) benchmark, t=',num2str(Ts3D(it))]);
%     title(['HeLaZ,$(P,J) =(',num2str(PMAXI),', ', num2str(JMAXI),'$), $\nu =',num2str(NU),...
%         '$, $\epsilon = ',num2str(eps),'$, $k_y = ', num2str(ky( iky)),'$, $q =',num2str(q0),'$, $s = ', num2str(shear),'$, $R_N = ', ...
%         num2str(K_N),'$, $R_{T_i} = ', num2str(K_T),'$, $N_z =',num2str(Nz),'$']);
    %set(gca,'Yscale','log')
    %
    
%     %Check symmetry of the mode at the outter mid plane
%     figure; hold all;
%     right = phib_real(midpoint+1:end);
%     left = fliplr(phib_real(1:midpoint-1)');
%     up_down_symm  = right(1:end) - left(1:end-1)';
%     %plot(b_angle(midpoint+1:end)/pi,phib_real(midpoint+1:end),'-xb');
%     plot(b_angle(midpoint+1:end)/pi,up_down_symm ,'-xb');
    %plot(abs(b_angle(1:midpoint-1))/pi,phib_real(1:midpoint-1),'-xb');
    %
    %
    % figure; hold all
    % plot(b_angle/pi, phib_imag.^2 + phib_real.^2 ,'xk');
    % %set(gca,'Yscale','log')
    % xlabel('$\chi / \pi$')
    % ylabel('$\phi_B (\chi)$');
    % title(['$(P,J) =(',num2str(pmax),', ', num2str(jmax),'$), $\nu =',num2str(nu),'$, $\epsilon = ',num2str(epsilon),'$, $q =',num2str(safety_fac),'$, $s = ', num2str(shear),'$, $k_y =',num2str(ky),'$']);
end
