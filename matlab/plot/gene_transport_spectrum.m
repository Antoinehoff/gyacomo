% folder = '/misc/gene_results/NL_Zpinch_Kn_1.8_eta_0.25_nuSG_5e-2_mu_1e-2_SGDK_36x20/';
% folder = '/misc/gene_results/HP_fig_2c_mu_5e-2/';
folder = '/misc/gene_results/HP_fig_2b_mu_5e-2/';
% folder = '/misc/gene_results/HP_fig_2c_gyroLES/';
% folder = '/misc/gene_results/NL_Zpinch_Kn_1.8_eta_0.25_nuSG_5e-2_mu_1e-2_SGDK_128x32x48x32/';

file = 'coord.dat.h5';
vp = h5read([folder,file],'/coord/vp');
mu = h5read([folder,file],'/coord/mu');
kx = h5read([folder,file],'/coord/kx');
ky = h5read([folder,file],'/coord/ky');
z  = h5read([folder,file],'/coord/z');

[KY,~] = meshgrid(ky,kx);

file = 'field.dat.h5';
time  = h5read([folder,file],'/field/time');

KDIM = 'ky';
TIMES = 2500:100:3500;
switch KDIM
    case 'kx'
        sdim  = 2;
        k     = kx;
        xname = '$k_x$';
        yname = '$\sum_{k_y}|\Gamma_k|$';
        shiftx = @(x) x(1:numel(kx)/2);
        shifty = @(x) x(1:numel(kx)/2);
    case 'ky'
        sdim  = 1;
        k     = ky;
        xname = '$k_y$';
        yname = '$\sum_{k_x}|\Gamma_k|$';
        shiftx = @(x) x;
        shifty = @(x) x(1:numel(ky));
end

colors = jet(numel(TIMES));
i_ = 1;
figure
for T = TIMES
[~, it] = min(abs(time-T));

v_x = h5read([folder,   'field.dat.h5'],[    '/field/phi/',sprintf('%10.10d',it-1)]);
v_x = v_x.real + 1i*v_x.imaginary; v_x = squeeze(v_x(:,:,1));
v_x = -1i*KY.*v_x;

n_i = h5read([folder,'mom_ions.dat.h5'],['/mom_ions/dens/',sprintf('%10.10d',it-1)]);
n_i = n_i.real + 1i*n_i.imaginary; n_i = squeeze(n_i(:,:,1));
% topclr_ = abs(n_i);

% v_x = computeFXYZ(v_x);
% n_i = computeFXYZ(n_i);
v_x = real(ifft2(v_x));
n_i = real(ifft2(n_i));

Gx  = n_i.*v_x;

Gx  = fft2(Gx);
% topclr_ = abs(Gx);
end
[~, it0] = min(abs(time-T(1)));
[~, it1] = min(abs(time-T(end)));
Gx = mean(Gx,3);

subplot(1,2,1)
    sdim  = 2;
    k     = kx;
    shiftx = @(x) x;%(1:numel(kx)/2);
    shifty = @(x) x;%(1:numel(kx)/2);
    Y  = abs(Gx(:,ky==0));
    Y  = Y/max(Y);
    semilogy(shiftx(k),shifty(Y),'DisplayName',['$t=',num2str(time(it)),'$'],...
            'Color',colors(i_,:)); hold on;

subplot(1,2,2)
    sdim  = 1;
    k     = ky;
    shiftx = @(x) x;
    shifty = @(x) x;%(1:numel(ky));
    Y  = abs(Gx(kx==0,:));
    Y  = Y/max(Y);
    semilogy(shiftx(k),shifty(Y),'DisplayName',['$t=',num2str(time(it)),'$'],...
            'Color',colors(i_,:)); hold on;
i_ = i_+1;

subplot(1,2,1)
    xname = '$k_x$';
    yname = '$\sum_{k_y}|\Gamma_k|$';
    title('Gene $k_x$ transport spectrum'); legend('show','Location','eastoutside');
    xlabel(xname); ylabel(yname)
    
subplot(1,2,2)
    xname = '$k_y$';
    yname = '$\sum_{k_x}|\Gamma_k|$';
    title('Gene $k_y$ transport spectrum'); legend('show','Location','eastoutside');
    xlabel(xname); ylabel(yname)
    
    figure; pclr = pcolor(topclr_); set(pclr,'EdgeColor','none');
