function [FIGURE, wkykx, ekykx] = mode_growth_meter(DATA,OPTIONS)

NORMALIZED = OPTIONS.NORMALIZED;
Nma   = OPTIONS.NMA; %Number moving average
d     = OPTIONS.fftz.flag;  % To add spectral evolution of z (useful for 3d zpinch)
switch OPTIONS.FIELD
    case 'Ni00'
        FIELD = DATA.Ni00;
        fname = 'N^{00}';
    case 'phi'
        FIELD = DATA.PHI;
        fname = '\phi';
end
if numel(size(FIELD)) == 3
        field = squeeze(FIELD);
        zstrcomp = 'z=0';
else
    switch OPTIONS.iz
        case 'avg'
            field = reshape(mean(FIELD,3),DATA.grids.Nky,DATA.grids.Nkx,numel(DATA.Ts3D));
            zstrcomp = 'z-avg';
        otherwise
            field = reshape(FIELD(:,:,OPTIONS.iz,:),DATA.grids.Nky,DATA.grids.Nkx,numel(DATA.Ts3D));
            zstrcomp = ['z=',num2str(DATA.grids.z(OPTIONS.iz))];
    end
end

FRAMES = zeros(size(OPTIONS.TIME));
for i = 1:numel(OPTIONS.TIME)
    [~,FRAMES(i)] =min(abs(OPTIONS.TIME(i)-DATA.Ts3D));
end
FRAMES = unique(FRAMES);
t  = DATA.Ts3D(FRAMES);

% time window where we measure the growth
TW = [OPTIONS.KY_TW; OPTIONS.KX_TW];

[~,ikzf] = max(mean(squeeze(abs(field(1,:,FRAMES))),2));

% Amplitude ratio method to measure growth rates and freq.
[~,it1] = min(abs(t-TW(1,1)));
[~,it2] = min(abs(t-TW(1,2)));
[wkykx, ekykx] = compute_growth_rates(FIELD(:,:,:,it1:it2),DATA.Ts3D(it1:it2));

FIGURE = struct();
if OPTIONS.SHOWFIG
FIGURE.fig = figure; %set(gcf, 'Position',  [100 100 1200 700])
FIGURE.FIGNAME = 'mode_growth_meter';
for i = 1:2
    MODES_SELECTOR = i; %(1:kx=0; 2:ky=0)

    [~,it1] = min(abs(t-TW(i,1)));
    [~,it2] = min(abs(t-TW(i,2)));
    
    if MODES_SELECTOR == 2
        switch OPTIONS.ik
            case 'sum'
                plt_ = @(f,ik) movmean(squeeze(sum(abs((f(:,ik,FRAMES))),1)),Nma);
                MODESTR = '$\sum_{k_y}$';
                omega= @(ik) wkykx(1,:);
                err  = @(ik) ekykx(1,:);
            case 'max'
                plt_ = @(f,ik) movmean(squeeze(max(abs((f(:,ik,FRAMES))),[],1)),Nma);
                MODESTR = '$\max_{k_y}$';
                omega= @(ik) wkykx(1,:);
                err  = @(ik) ekykx(1,:);
            otherwise
                plt_ = @(f,ik) movmean(abs(squeeze(f(OPTIONS.ik,ik,FRAMES))),Nma);
                MODESTR = ['$k_y=$',num2str(DATA.grids.ky(OPTIONS.ik))];
                omega= @(ik) wkykx(OPTIONS.ik,:);
                err  = @(ik) ekykx(OPTIONS.ik,:);
        end
        kstr = 'k_x';
        % Number max of modes to plot is kx>0 (1/2) of the non filtered modes (2/3)
        Nmax = ceil(DATA.grids.Nkx*1/3);
        k = DATA.grids.kx;
    elseif MODES_SELECTOR == 1
        switch OPTIONS.ik
            case 'sum'
                plt_ = @(x,ik) movmean(squeeze(sum(abs((x(ik,:,FRAMES))),2)),Nma);
                MODESTR = '$\sum_{k_x}$';
                omega= @(ik) wkykx(:,1);
                err  = @(ik) ekykx(:,1);
            case 'max'
                plt_ = @(x,ik) movmean(squeeze(max(abs((x(ik,:,FRAMES))),[],2)),Nma);
                MODESTR = '$\max_{k_x}$';
                omega= @(ik) wkykx(:,1);
                err  = @(ik) ekykx(:,1);
            otherwise
                plt_ = @(x,ik) movmean(abs(squeeze(x(ik,OPTIONS.ik,FRAMES))),Nma);
                MODESTR = ['$k_x=$',num2str(DATA.grids.kx(OPTIONS.ik))];
                omega= @(ik) wkykx(:,OPTIONS.ik);
                err  = @(ik) ekykx(:,OPTIONS.ik);
       end
        kstr = 'k_y';
        % Number max of modes to plot is ky>0 (1/1) of the non filtered modes (2/3)
        Nmax = ceil(DATA.grids.Nky*2/3);
        k = DATA.grids.ky;
    end 
    if NORMALIZED
        plt = @(x,ik) plt_(x,ik)./max(plt_(x,ik));
    else
        plt = @(x,ik) plt_(x,ik);
    end

    MODES = 1:numel(k);

    gamma = MODES;
    amp   = MODES;
    i_=1;
    % Alternative method to measure growth
    for ik = MODES
        to_measure = log(plt(field,ik));
        gr = polyfit(t(it1:it2),to_measure(it1:it2),1);
        gamma(i_) = gr(1);
        amp(i_)   = gr(2);
        i_=i_+1;
    end
    mod2plot = 2:min(Nmax,OPTIONS.NMODES+1);
    clr_ = jet(numel(mod2plot));
    %plot
%     subplot(2+d,3,1+3*(i-1))
    FIGURE.axes(1+3*(i-1)) = subplot(2+d,3,1+3*(i-1),'parent',FIGURE.fig);
        [YY,XX] = meshgrid(t,fftshift(k));
        pclr = pcolor(XX,YY,abs(plt(fftshift(field,MODES_SELECTOR),1:numel(k))));set(pclr, 'edgecolor','none');  hold on;
        set(gca,'YDir','normal')
    %     xlim([t(1) t(end)]); %ylim([1e-5 1])
        xlabel(['$',kstr,'\rho_s$']); ylabel('$t c_s /R_0$');
        title([MODESTR,', ',zstrcomp])  

%     subplot(2+d,3,2+3*(i-1))
    FIGURE.axes(2+3*(i-1)) = subplot(2+d,3,2+3*(i-1),'parent',FIGURE.fig);
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
        xlabel('$t c_s /R_0$'); ylabel(['$|',fname,'_{',kstr,'}|$']);
        title('Measure time window')

%     subplot(2+d,3,3+3*(i-1))
    FIGURE.axes(3+3*(i-1)) = subplot(2+d,3,3+3*(i-1),'parent',FIGURE.fig);
    % yyaxis("left")
    if OPTIONS.GOK2
        x_ = k(MODES);
        y_ = real(omega(ik))./k(MODES).^2;
        e_ = real(err(ik))./k(MODES).^2;
        lb_= '\gamma/k^2';
    else
        x_ = k(MODES);
        y_ = real(omega(ik));
        e_ = real(err(ik));
        lb_= '\gamma';
    end
        errorbar(x_,y_,e_,'-ok',...
                'LineWidth',1.5,'MarkerSize',10,...
                'DisplayName',...
                ['$\gamma$, (',num2str(DATA.inputs.PMAX),',',num2str(DATA.inputs.JMAX),')']); 
        hold on;
        for i_ = 1:numel(mod2plot)
            plot(k(MODES(mod2plot(i_))),gamma(mod2plot(i_)),'x','color',clr_(i_,:));
        end
        if MODES_SELECTOR == 2
            plot(k(ikzf),gamma(ikzf),'ok');
        end
        % ylabel('$\gamma$');
    % yyaxis("right")
        errorbar(k(MODES),imag(omega(ik)),imag(err(ik)),'*k',...
                'LineWidth',1.5,'MarkerSize',10,...
                'DisplayName',...
                ['$\omega$, (',num2str(DATA.inputs.PMAX),',',num2str(DATA.inputs.JMAX),')']); 
        hold on;
        ylabel(['$',lb_,',\omega$']);
        xlabel(['$',kstr,'\rho_s$']);
        title('Growth rates')
end

if d
    [~,ikx] = min(abs(DATA.grids.kx-OPTIONS.fftz.kx));
    [~,iky] = min(abs(DATA.grids.ky-OPTIONS.fftz.ky));
    sz_=size(FIELD);nkz = sz_(3)/2;
    k = [(0:nkz/2), (-nkz/2+1):-1]/DATA.grids.Npol;
    % Spectral treatment of z-dimension
    Y = fft(FIELD(iky,ikx,:,:),[],3);
    phi = squeeze(Y(1,1,2:2:end,:)); 
    
    gamma = zeros(nkz);
    for ikz = 1:nkz
        to_measure = squeeze(log(abs(phi(ikz,it1:it2))));
        tmp = polyfit(t(it1:it2),to_measure(:),1);
        if ~(isnan(tmp(1)) || isinf(tmp(1)))
            gamma(ikz) = tmp(1);
        end
    end
     %plot
    MODES    = 1:numel(k);
    mod2plot = [2:2:min(nkz,OPTIONS.NMODES+1)];
    clr_     = jet(numel(mod2plot));
    subplot(3,3,7)
        [YY,XX] = meshgrid(t,fftshift(k));
        pclr = pcolor(XX,YY,abs(phi));set(pclr, 'edgecolor','none');  hold on;
        set(gca,'YDir','normal')
    %     xlim([t(1) t(end)]); %ylim([1e-5 1])
        xlabel('$k_\parallel$'); ylabel('$t c_s /\rho_s$');
        title(['$k_x=$',num2str(DATA.grids.kx(ikx)),', $k_y=$',num2str(DATA.grids.ky(iky))]);  

    subplot(3,3,8)
        for i_ = 1:numel(mod2plot)
            im_ = MODES(mod2plot(i_));
            semilogy(t,abs(phi(im_,:)),'color',clr_(i_,:)); hold on;
        end
        %plot the time window where the gr are measured
        plot(t(it1)*[1 1],[1e-10 1],'--k')
        plot(t(it2)*[1 1],[1e-10 1],'--k')
        xlim([t(1) t(end)]); %ylim([1e-5 1])
        xlabel('$t c_s /\rho_s$'); ylabel(['$|\phi_{',kstr,'}|$']);
        title('Measure time window')

    subplot(3,3,9)
        plot(k(MODES),gamma,...
                'DisplayName',['(',num2str(DATA.inputs.PMAX),',',num2str(DATA.inputs.JMAX),')']); hold on;
        for i_ = 1:numel(mod2plot)
            plot(k(MODES(mod2plot(i_))),gamma(mod2plot(i_)),'x','color',clr_(i_,:));
        end
        if MODES_SELECTOR == 2
            plot(k(ikzf),gamma(ikzf),'ok');
        end
        xlabel(['$k_\parallel$']); ylabel('$\gamma$');
        title('Growth rates')   
    
end
top_title(DATA.paramshort)
end
end