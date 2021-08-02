
%%
if 0
%% Hermite energy spectra
% tf = Ts2D(end-3);
skip = 1;
fig = figure; FIGNAME = ['hermite_spectrum_',PARAMS];set(gcf, 'Position',  [100, 100, 1800, 600]);
plt = @(x) squeeze(x);
N10 = floor(Nji/10);
Nrow = (1+N10);
Ncol = ceil(Nji/Nrow);
for ij = 1:Nji
    subplot(Nrow,Ncol,ij)
    for it5 = 1:skip:Ns5D
        alpha = it5*1.0/Ns5D;
        loglog(Pi(1:2:end),plt(epsilon_i_pj(1:2:end,ij,it5)),...
            'color',(1-alpha)*[0.8500, 0.3250, 0.0980]+alpha*[0, 0.4470, 0.7410],...
            'DisplayName',['t=',num2str(Ts5D(it5))]); hold on;
    end
    grid on; ylim([1e0,1e10]);
    xlabel('$p$');
    TITLE = ['$\sum_{kx,ky} |N_i^{p',num2str(Ji(ij)),'}|^2$']; title(TITLE);
end
save_figure
end
%%
if 0
%% Laguerre energy spectra
% tf = Ts2D(end-3);
skip = 1;
fig = figure; FIGNAME = ['laguerre_spectrum_',PARAMS];set(gcf, 'Position',  [100, 100, 1800, 600]);
plt = @(x) squeeze(x);
N10 = floor(Npi/10);
Nrow = ceil((1+N10)/2);
Ncol = ceil(Npi/Nrow/2);
for ip = 1:2:Npi
    subplot(Nrow,Ncol,ip/2+1)
    for it5 = 1:skip:Ns5D
        alpha = it5*1.0/Ns5D;
        loglog(Ji,plt(epsilon_i_pj(ip,:,it5)),...
            'color',(1-alpha)*[0.8500, 0.3250, 0.0980]+alpha*[0, 0.4470, 0.7410],...
            'DisplayName',['t=',num2str(Ts5D(it5))]); hold on;
        grid on;
        xlabel('$j$'); ylim([1e-20,1e10]);
        TITLE = ['$\sum_{kx,ky} |N_i^{',num2str(Pi(ip)),'j}|^2$']; title(TITLE);
    end
end
save_figure
end

%%
no_AA     = (2:floor(2*Nkx/3));
tKHI      = 100;
[~,itKHI] = min(abs(Ts2D-tKHI));
after_KHI = (itKHI:Ns2D);
if 0
%% Phi frequency space time diagram at ky=ky(iky)
ky_ = 0.0;
[~,iky] = min(abs(ky-ky_));
fig = figure; FIGNAME = ['phi_freq_diag_',PARAMS];set(gcf, 'Position',  [100, 100, 500, 400]);
        [TY,TX] = meshgrid(Ts2D(after_KHI),kx(no_AA));
        pclr = pcolor(TX,TY,(squeeze(abs(PHI(no_AA,iky,(after_KHI)))))); set(pclr, 'edgecolor','none'); colorbar;
        caxis([0,10000]); colormap hot
        ylabel('$t c_s/R$'), xlabel('$0<k_r<2/3 k_r^{\max}$')
        legend(['$|\tilde\phi(k_z=',num2str(ky_),')|$'])
        title('Spectrogram of $\phi$')
end