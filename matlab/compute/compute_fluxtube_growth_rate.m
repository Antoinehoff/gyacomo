function [ linear_gr ] = compute_fluxtube_growth_rate(DATA, TRANGE, PLOT)

% Remove AA part
if DATA.Nx > 1
    ikxnz = abs(DATA.PHI(1,:,1,1)) > 0;
else
    ikxnz = abs(DATA.PHI(1,:,1,1)) > -1;
end
ikynz = abs(DATA.PHI(:,1,1,1)) > 0;

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
if PLOT >0
       figure
       plot(linear_gr.ky,linear_gr.g_ky(:,end),'-o','DisplayName','$Re(\omega_{k_y})$'); hold on;
       plot(linear_gr.ky,linear_gr.w_ky(:,end),'-o','DisplayName','$Im(\omega_{k_y})$'); hold on;
       plot(linear_gr.ky,linear_gr.ce  (:,end),'-o','DisplayName','$\epsilon$'); hold on;
       xlim([min(linear_gr.ky) max(linear_gr.ky)]);
       xlabel('$k_y$');
       legend('show');
       title(DATA.param_title);
       if PLOT > 1
           [YY,XX] = meshgrid(kys,t(its+1:ite));
           figure
              subplot(311)
%            imagesc(t(its:ite),kys,real(w_ky)); set(gca,'YDir','normal'); 
           pclr = pcolor(XX',YY',real(w_ky));    set(pclr, 'edgecolor','none'); 
           xlabel('$t$'); ylabel('$k_y$');
           title('Re($\omega$)')
           
              subplot(312)
           pclr = pcolor(XX',YY',imag(w_ky));    set(pclr, 'edgecolor','none'); 
           xlabel('$t$'); ylabel('$k_y$');
           title('Im($\omega$)')
           
              subplot(313)
           pclr = pcolor(XX',YY',abs(w_ky));    set(pclr, 'edgecolor','none'); 
           xlabel('$t$'); ylabel('$k_y$');
           title('|$\omega$|')
       end

end

end