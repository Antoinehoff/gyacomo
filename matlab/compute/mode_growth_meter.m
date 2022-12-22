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

% time window where we measure the growth
it1 = floor(numel(t)/2);
it2 = numel(t);

[~,ikzf] = max(mean(squeeze(abs(field(1,:,FRAMES))),2));

FIGURE.fig = figure; set(gcf, 'Position',  [100 100 1200 700])
FIGURE.FIGNAME = 'mode_growth_meter';
for i = 1:2
    MODES_SELECTOR = i; %(2:Zonal, 1: NZonal)

    if MODES_SELECTOR == 2
        switch OPTIONS.ik
            case 'sum'
                plt_ = @(x,ik) movmean(squeeze(sum(abs((x(:,ik,FRAMES))),1)),Nma);
                MODESTR = '$\sum_{k_y}$';
            case 'max'
                plt_ = @(x,ik) movmean(squeeze(max(abs((x(:,ik,FRAMES))),[],1)),Nma);
                MODESTR = '$\max_{k_y}$';
            otherwise
                plt_ = @(x,ik) movmean(abs(squeeze(x(OPTIONS.ik,ik,FRAMES))),Nma);
                MODESTR = ['$k_y=$',num2str(DATA.ky(OPTIONS.ik))];
        end
        kstr = 'k_x';
        % Number max of modes to plot is kx>0 (1/2) of the non filtered modes (2/3)
        Nmax = ceil(DATA.Nkx*1/3);
        k = DATA.kx;
    elseif MODES_SELECTOR == 1
        switch OPTIONS.ik
            case 'sum'
                plt_ = @(x,ik) movmean(squeeze(sum(abs((x(ik,:,FRAMES))),2)),Nma);
                MODESTR = '$\sum_{k_x}$';
            case 'max'
                plt_ = @(x,ik) movmean(squeeze(max(abs((x(ik,:,FRAMES))),[],2)),Nma);
                MODESTR = '$\max_{k_x}$';
            otherwise
                plt_ = @(x,ik) movmean(abs(squeeze(x(ik,OPTIONS.ik,FRAMES))),Nma);
                MODESTR = ['$k_x=$',num2str(DATA.kx(OPTIONS.ik))];
        end
        kstr = 'k_y';
        % Number max of modes to plot is ky>0 (1/1) of the non filtered modes (2/3)
        Nmax = ceil(DATA.Nky*2/3);
        k = DATA.ky;
    end 
    if NORMALIZED
        plt = @(x,ik) plt_(x,ik)./max(plt_(x,ik));
    else
        plt = @(x,ik) plt_(x,ik);
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
        to_measure = log(plt(field,ik));
        gr = polyfit(t(it1:it2),to_measure(it1:it2),1);
        gamma(i_) = gr(1);
        amp(i_)   = gr(2);
        i_=i_+1;
    end
    mod2plot = [2:min(Nmax,OPTIONS.NMODES+1)];
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
        %plot the time window where the gr are measured
        plot(t(it1)*[1 1],[1e-10 1],'--k')
        plot(t(it2)*[1 1],[1e-10 1],'--k')
        xlim([t(1) t(end)]); %ylim([1e-5 1])
        xlabel('$t c_s /\rho_s$'); ylabel(['$|\phi_{',kstr,'}|$']);
        title('Measure time window')

    subplot(2,3,3+3*(i-1))
        plot(k(MODES),gamma,...
                'DisplayName',['(',num2str(DATA.Pmaxi-1),',',num2str(DATA.Jmaxi-1),')']); hold on;
        for i_ = 1:numel(mod2plot)
            plot(k(MODES(mod2plot(i_))),gamma(mod2plot(i_)),'x','color',clr_(i_,:));
        end
        if MODES_SELECTOR == 2
            plot(k(ikzf),gamma(ikzf),'ok');
        end
        xlabel(['$',kstr,'\rho_s$']); ylabel('$\gamma$');
        title('Growth rates')
end
end