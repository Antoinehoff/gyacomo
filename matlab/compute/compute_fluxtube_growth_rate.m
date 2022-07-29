function [ linear_gr ] = compute_fluxtube_growth_rate(DATA, TRANGE, PLOT)

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

[~,its] = min(abs(t-TRANGE(1)));
[~,ite] = min(abs(t-TRANGE(end)));

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

linear_gr.trange = t(its:ite);
linear_gr.g_ky   = real(w_ky(Is,:));
linear_gr.w_ky   = imag(w_ky(Is,:));
linear_gr.ce     = abs(ce);
linear_gr.ky     = kys;
linear_gr.std_g = std (real(w_ky(Is,:)),[],2);
linear_gr.avg_g = mean(real(w_ky(Is,:)),2);
linear_gr.std_w = std (imag(w_ky(Is,:)),[],2);
linear_gr.avg_w = mean(imag(w_ky(Is,:)),2);

if PLOT >0
   figure
if PLOT > 1
    subplot(1,2,1)
end

       plot(linear_gr.ky,linear_gr.g_ky(:,end),'-o','DisplayName','$Re(\omega_{k_y})$'); hold on;
       plot(linear_gr.ky,linear_gr.w_ky(:,end),'--*','DisplayName','$Im(\omega_{k_y})$'); hold on;
       plot(linear_gr.ky,linear_gr.ce  (:,end),'x','DisplayName','$\epsilon$'); hold on;

       errorbar(linear_gr.ky,linear_gr.avg_g,linear_gr.std_g,...
           '-o','DisplayName','$Re(\omega_{k_y})$',...
           'LineWidth',1.5); hold on;
       errorbar(linear_gr.ky,linear_gr.avg_w,linear_gr.std_w,...
           '--*','DisplayName','$Im(\omega_{k_y})$',...
           'LineWidth',1.5); hold on;

%        xlim([min(linear_gr.ky) max(linear_gr.ky)]);
       xlabel('$k_y$');
       legend('show');
       title(DATA.param_title);
       
if PLOT > 1
    if PLOT == 2
    subplot(1,2,2)
    elseif PLOT == 3
        subplot(2,2,2)
    end
    plot(DATA.Ts3D(its+1:ite),linear_gr.g_ky(Is,:)); hold on;
    plot(DATA.Ts3D(its+1:ite),linear_gr.w_ky(Is,:));
    xlabel('t'); ylabel('$\gamma(t),\omega(t)$'); xlim([DATA.Ts3D(1) DATA.Ts3D(end)]);
end

if PLOT > 2
    xlabel([]); xticks([]);
    subplot(2,2,4)
    semilogy(DATA.Ts3D,squeeze(abs(DATA.PHI(2,1,DATA.Nz/2,:)))); hold on;
    xlabel('t'); ylabel('$|\phi_{ky}|(t)$')

end

end