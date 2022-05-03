function [FIGURE] = mode_growth_meter(DATA,OPTIONS)

NORMALIZED = OPTIONS.NORMALIZED;
Nma   = OPTIONS.NMA; %Number moving average
t  = OPTIONS.TIME;
iz = OPTIONS.iz;
[~,ikzf] = max(squeeze(mean(abs(squeeze(DATA.PHI(1,:,1,:))),2)));

FRAMES = zeros(size(OPTIONS.TIME));

for i = 1:numel(OPTIONS.TIME)
    [~,FRAMES(i)] =min(abs(OPTIONS.TIME(i)-DATA.Ts3D));
end

FIGURE.fig = figure; set(gcf, 'Position',  [100 100 1200 700])
FIGURE.FIGNAME = 'mode_growth_meter';
for i = 1:2
MODES_SELECTOR = i; %(1:Zonal, 2: NZonal, 3: ky=kx)

if MODES_SELECTOR == 1
    if NORMALIZED
        plt = @(x,ik) movmean(abs(squeeze(x(1,ik,iz,FRAMES)))./max(abs(squeeze(x(1,ik,iz,FRAMES)))),Nma);
    else
        plt = @(x,ik) movmean(abs(squeeze(x(1,ik,iz,FRAMES))),Nma);
    end
    kstr = 'k_x';
    k = DATA.kx;
    MODESTR = 'Zonal modes';
elseif MODES_SELECTOR == 2
    if NORMALIZED
        plt = @(x,ik) movmean(abs(squeeze(x(ik,1,iz,FRAMES)))./max(abs(squeeze(x(ik,1,iz,FRAMES)))),Nma);
    else
        plt = @(x,ik) movmean(abs(squeeze(x(ik,1,iz,FRAMES))),Nma);
    end
    kstr = 'k_y';
    k = DATA.ky;
    MODESTR = 'NZ modes';
elseif MODES_SELECTOR == 3
    if NORMALIZED
        plt = @(x,ik) movmean(abs(squeeze(x(ik,ik,iz,FRAMES)))./max(abs(squeeze(x(ik,ik,iz,FRAMES)))),Nma);
    else
        plt = @(x,ik) movmean(abs(squeeze(x(ik,ik,iz,FRAMES))),Nma);
    end
    kstr = 'k_y=k_x';
    k = DATA.ky;
    MODESTR = 'Diag modes';
end 

MODES = 1:numel(k);
% MODES = zeros(size(OPTIONS.K2PLOT));
% for i = 1:numel(OPTIONS.K2PLOT)
%     [~,MODES(i)] =min(abs(OPTIONS.K2PLOT(i)-k));
% end


% plt = @(x,ik) abs(squeeze(x(1,ik,iz,FRAMES)))./max(abs(squeeze(x(1,ik,iz,FRAMES))));

gamma = MODES;
amp   = MODES;
i_=1;
for ik = MODES
    gr = polyfit(t',log(plt(DATA.PHI,ik)),1);
    gamma(i_) = gr(1);
    amp(i_)   = gr(2);
    i_=i_+1;
end

%plot
subplot(2,3,1+3*(i-1))
    [YY,XX] = meshgrid(t,fftshift(k,1));
    pclr = pcolor(XX,YY,abs(plt(fftshift(DATA.PHI,MODES_SELECTOR),1:numel(k))));set(pclr, 'edgecolor','none');  hold on;
    set(gca,'YDir','normal')
%     xlim([t(1) t(end)]); %ylim([1e-5 1])
    xlabel(['$',kstr,'\rho_s$']); ylabel('$t c_s /\rho_s$');
    title(MODESTR)  
    
subplot(2,3,2+3*(i-1))
    mod2plot = [2:OPTIONS.NMODES+1];
    for i_ = mod2plot
        semilogy(t,plt(DATA.PHI,MODES(i_))); hold on;
%         semilogy(t,exp(gamma(i_).*t+amp(i_)),'--k')
    end
    if MODES_SELECTOR == 1
        semilogy(t,plt(DATA.PHI,ikzf),'--k');
    end
    xlim([t(1) t(end)]); %ylim([1e-5 1])
    xlabel('$t c_s /\rho_s$'); ylabel(['$|\phi_{',kstr,'}|$']);
    title('Measure time window')
    
subplot(2,3,3+3*(i-1))
    plot(k(MODES),gamma); hold on;
    plot(k(MODES(mod2plot)),gamma(mod2plot),'x')
    if MODES_SELECTOR == 1
        plot(k(ikzf),gamma(ikzf),'ok');
    end
    xlabel(['$',kstr,'\rho_s$']); ylabel('$\gamma$');
    title('Growth rates')
end
end