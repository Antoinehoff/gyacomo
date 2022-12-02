function [ linear_gr ] = compute_fluxtube_growth_rate(DATA, OPTIONS)

% Remove AA part
if DATA.Nx > 1
    ikxnz = abs(DATA.PHI(1,:,1,1)) > 0;
else
    ikxnz = abs(DATA.PHI(1,:,1,1)) > -1;
end
ikynz = (abs(DATA.PHI(:,1,1,1)) > 0);
ikynz = logical(ikynz .* (DATA.ky > 0));
phi = DATA.PHI(ikynz,ikxnz,:,:);
t   = DATA.Ts3D;

[~,its] = min(abs(t-OPTIONS.TRANGE(1)));
[~,ite] = min(abs(t-OPTIONS.TRANGE(end)));

w_ky = zeros(sum(ikynz),ite-its);
ce   = zeros(sum(ikynz),ite-its);

is = 1;
for it = its+1:ite
    phi_n   = phi(:,:,:,it); 
    phi_nm1 = phi(:,:,:,it-1);
    dt      = t(it)-t(it-1);
    ZS      = sum(sum(phi_nm1,2),3);
   
    wl         = log(phi_n./phi_nm1)/dt;
    w_ky(:,is) = squeeze(sum(sum(wl.*phi_nm1,2),3)./ZS);
    
    for iky = 1:numel(w_ky(:,is))
        ce(iky,is)   = abs(sum(sum(abs(w_ky(iky,is)-wl(iky,:,:)).^2.*phi_nm1(iky,:,:),2),3)./ZS(iky,:,:));
    end
    is = is + 1;
end


[kys, Is] = sort(DATA.ky(ikynz));

linear_gr.OPTIONS.TRANGE = t(its:ite);
linear_gr.g_ky   = real(w_ky(Is,:));
linear_gr.w_ky   = imag(w_ky(Is,:));
linear_gr.ce     = abs(ce);
linear_gr.ky     = kys;
linear_gr.std_g = std (real(w_ky(Is,:)),[],2);
linear_gr.avg_g = mean(real(w_ky(Is,:)),2);
linear_gr.std_w = std (imag(w_ky(Is,:)),[],2);
linear_gr.avg_w = mean(imag(w_ky(Is,:)),2);

if OPTIONS.NPLOTS >0
   figure
if OPTIONS.NPLOTS > 1
    subplot(1,2,1)
end
    x_ = linear_gr.ky;
    plt = @(x) x;
    OVERK = '';
    if OPTIONS.GOK == 1
        plt = @(x) x./x_;
        OVERK = '/$k_y$';
    elseif OPTIONS.GOK == 2
        plt = @(x) x.^2./x_.^3;
        OVERK = '/$k_y$';
    end
       plot(x_,plt(linear_gr.g_ky(:,end)),'-o','DisplayName',['$Re(\omega_{k_y})$',OVERK]); hold on;
       plot(x_,plt(linear_gr.w_ky(:,end)),'--*','DisplayName',['$Im(\omega_{k_y})$',OVERK]); hold on;
       plot(x_,plt(linear_gr.ce  (:,end)),'x','DisplayName',['$\epsilon$',OVERK]); hold on;

       errorbar(linear_gr.ky,plt(linear_gr.avg_g),plt(linear_gr.std_g),...
           '-o','DisplayName','$Re(\omega_{k_y})$',...
           'LineWidth',1.5); hold on;
       errorbar(linear_gr.ky,plt(linear_gr.avg_w),plt(linear_gr.std_w),...
           '--*','DisplayName','$Im(\omega_{k_y})$',...
           'LineWidth',1.5); hold on;

%        xlim([min(linear_gr.ky) max(linear_gr.ky)]);
       xlabel('$k_y$');
       legend('show');
       title(DATA.param_title);
       
if OPTIONS.NPLOTS > 1
    if OPTIONS.NPLOTS == 2
    subplot(1,2,2)
    elseif OPTIONS.NPLOTS == 3
        subplot(2,2,2)
    end
    plot(DATA.Ts3D(its+1:ite),linear_gr.g_ky(Is,:)); hold on;
    plot(DATA.Ts3D(its+1:ite),linear_gr.w_ky(Is,:));
    xlabel('t'); ylabel('$\gamma(t),\omega(t)$'); xlim([DATA.Ts3D(1) DATA.Ts3D(end)]);
end

if OPTIONS.NPLOTS > 2
    xlabel([]); xticks([]);
    subplot(2,2,4)
    [~,ikx0] = min(abs(DATA.kx));
    for i_ = 1:(DATA.Nkx+1)/2
        iky = 1 + (DATA.ky(1) == 0);
        ikx = ikx0 + (i_-1);
        semilogy(DATA.Ts3D,squeeze(abs(DATA.PHI(iky,ikx,DATA.Nz/2,:))),...
            'DisplayName',['$k_x=',num2str(DATA.kx(ikx)),'$, $k_y=',num2str(DATA.ky(iky)),'$']); 
        hold on;
        xlabel('t,'); ylabel('$|\phi_{ky}|(t)$')
    end
legend('show')
end

end