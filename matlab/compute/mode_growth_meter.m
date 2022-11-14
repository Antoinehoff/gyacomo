function [FIGURE] = mode_growth_meter(DATA,OPTIONS)

NORMALIZED = OPTIONS.NORMALIZED;
Nma   = OPTIONS.NMA; %Number moving average

switch OPTIONS.iz
    case 'avg'
        field = squeeze(mean(DATA.PHI,3));
    otherwise
        field = squeeze(DATA.PHI(:,:,OPTIONS.iz,:));
end

FRAMES = zeros(size(OPTIONS.TIME));
for i = 1:numel(OPTIONS.TIME)
    [~,FRAMES(i)] =min(abs(OPTIONS.TIME(i)-DATA.Ts3D));
end
FRAMES = unique(FRAMES);
t  = DATA.Ts3D(FRAMES);

[~,ikzf] = max(mean(squeeze(abs(field(1,:,FRAMES))),2));

FIGURE.fig = figure; set(gcf, 'Position',  [100 100 1200 700])
FIGURE.FIGNAME = 'mode_growth_meter';
for i = 1:2
    MODES_SELECTOR = i; %(2:Zonal, 1: NZonal, 3: ky=kx)

    if MODES_SELECTOR == 2
        if NORMALIZED
            plt = @(x,ik) movmean(abs(squeeze(x(1,ik,FRAMES)))./max(abs(squeeze(x(1,ik,FRAMES)))),Nma);
        else
            plt = @(x,ik) movmean(abs(squeeze(x(1,ik,FRAMES))),Nma);
        end
        kstr = 'k_x';
        k = DATA.kx;
        MODESTR = 'Zonal modes';
    elseif MODES_SELECTOR == 1
        if NORMALIZED
            plt = @(x,ik) movmean(abs(squeeze(x(ik,1,FRAMES)))./max(abs(squeeze(x(ik,1,FRAMES)))),Nma);
        else
            plt = @(x,ik) movmean(abs(squeeze(x(ik,1,FRAMES))),Nma);
        end
        kstr = 'k_y';
        k = DATA.ky;
        MODESTR = 'NZ modes';
    elseif MODES_SELECTOR == 3
        if NORMALIZED
            plt = @(x,ik) movmean(abs(squeeze(x(ik,ik,FRAMES)))./max(abs(squeeze(x(ik,ik,FRAMES)))),Nma);
        else
            plt = @(x,ik) movmean(abs(squeeze(x(ik,ik,FRAMES))),Nma);
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
        gr = polyfit(t,log(plt(field,ik)),1);
        gamma(i_) = gr(1);
        amp(i_)   = gr(2);
        i_=i_+1;
    end
    mod2plot = [2:OPTIONS.NMODES+1];
    clr_ = jet(numel(mod2plot));
    %plot
    subplot(2,3,1+3*(i-1))
        [YY,XX] = meshgrid(t,fftshift(k));
        pclr = pcolor(XX,YY,abs(plt(fftshift(field,MODES_SELECTOR),1:numel(k))));set(pclr, 'edgecolor','none');  hold on;
        set(gca,'YDir','normal')
    %     xlim([t(1) t(end)]); %ylim([1e-5 1])
        xlabel(['$',kstr,'\rho_s$']); ylabel('$t c_s /\rho_s$');
        title(MODESTR)  

    subplot(2,3,2+3*(i-1))
        for i_ = 1:numel(mod2plot)
            semilogy(t,plt(field,MODES(mod2plot(i_))),'color',clr_(i_,:)); hold on;
    %         semilogy(t,exp(gamma(i_).*t+amp(i_)),'--k')
        end
        if MODES_SELECTOR == 2
            semilogy(t,plt(field,ikzf),'--k');
        end
        xlim([t(1) t(end)]); %ylim([1e-5 1])
        xlabel('$t c_s /\rho_s$'); ylabel(['$|\phi_{',kstr,'}|$']);
        title('Measure time window')

    subplot(2,3,3+3*(i-1))
        plot(k(MODES),gamma); hold on;
        for i_ = 1:numel(mod2plot)
            plot(k(MODES(mod2plot(i_))),gamma(mod2plot(i_)),'x','color',clr_(i_,:))
        end
        if MODES_SELECTOR == 2
            plot(k(ikzf),gamma(ikzf),'ok');
        end
        xlabel(['$',kstr,'\rho_s$']); ylabel('$\gamma$');
        title('Growth rates')
end
end