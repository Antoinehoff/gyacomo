function [ FIGURE ] = show_moments_spectrum( DATA, OPTIONS )


Pi = DATA.Pi;
Ji = DATA.Ji;
Nipj = sum(sum(sum(abs(DATA.Nipj),3),4),5);
Nipj = squeeze(Nipj);

if DATA.KIN_E
Pe = DATA.Pe;
Je = DATA.Je;
Nepj = sum(sum(sum(abs(DATA.Nepj),3),4),5);
Nepj = squeeze(Nepj);
end

FIGURE.fig = figure; FIGURE.FIGNAME = ['mom_spectrum_',DATA.PARAMS];
set(gcf, 'Position',  [100 50 1000 400])

if ~OPTIONS.ST
    t0 = OPTIONS.TIME(1);
    t1 = OPTIONS.TIME(end);
    [~,it0] = min(abs(t0-DATA.Ts5D));
    [~,it1] = min(abs(t1-DATA.Ts5D));
    Nipj = mean(Nipj(:,:,it0:it1),3);
    if DATA.K_E
    Nepj = mean(Nepj(:,:,it0:it1),3);
    end
    if numel(OPTIONS.TIME) == 1
        TITLE=['$t=',num2str(OPTIONS.TIME),'$'];
    else
        TITLE=['$t\in[',num2str(t0),',',num2str(t1),']$'];
    end 
    Nipj = squeeze(Nipj);

    ymini = min(Nipj); ymaxi = max(Nipj);

    if DATA.K_E
    Nepj = squeeze(Nepj);
    ymine = min(Nepj); ymaxe = max(Nepj);
    ymax  = max([ymaxi ymaxe]);
    ymin  = min([ymini ymine]);
    end
    if DATA.K_E
    subplot(121)
    end
    if ~OPTIONS.P2J
        for ij = 1:numel(Ji)
            name = ['$j=',num2str(Ji(ij)),'$'];
            semilogy(Pi,Nipj(:,ij),'o-','DisplayName',name); hold on;
        end
        xlabel('$p$'); 
    else
        for ij = 1:numel(Ji)
            name = ['$j=',num2str(Ji(ij)),'$'];
            semilogy(Pi+2*Ji(ij),Nipj(:,ij),'o-','DisplayName',name); hold on;
        end
        xlabel('$p+2j$')
    end
    ylabel(['$\sum_{kx,ky}|N_i^{pj}|$']);
    legend('show');
    title([TITLE,' He-La ion spectrum']);
    
    if DATA.K_E
    subplot(122)
    if ~OPTIONS.P2J
        for ij = 1:numel(Je)
            name = ['$j=',num2str(Je(ij)),'$'];
            semilogy(Pe,Nepj(:,ij),'o-','DisplayName',name); hold on;
        end
        xlabel('p'); 
    else
    for ij = 1:numel(Je)
        name = ['$j=',num2str(Je(ij)),'$'];
        semilogy(Pe+2*Je(ij),Nepj(:,ij),'o-','DisplayName',name); hold on;
    end
    xlabel('$p+2j$')
    end
    ylabel(['$\sum_{kx,ky}|N_e^{pj}|$']);
    legend('show');
    title([TITLE,' He-La elec. spectrum']);
    end
else
    [JJ,PP] = meshgrid(Ji,Pi);
    P2Ji = PP + 2*JJ;
    % form an axis of velocity ordered moments
    p2ji = unique(P2Ji);
    % weights to average
    weights = zeros(size(p2ji));
    % space time of moments amplitude wrt to p+2j degree
    Ni_ST = zeros(numel(p2ji),numel(DATA.Ts5D));
    % fill the st diagramm by averaging every p+2j=n moments together
    for ip = 1:numel(Pi)
        for ij = 1:numel(Ji)
            [~,ip2j] = min(abs(p2ji-(Pi(ip)+2*Ji(ij))));
            Ni_ST(ip2j,:) = Ni_ST(ip2j,:) + transpose(squeeze(Nipj(ip,ij,:)));
            weights(ip2j) = weights(ip2j) + 1;
        end
    end
    % doing the average
    for ip2j = 1:numel(p2ji)
        Ni_ST(ip2j,:) = Ni_ST(ip2j,:)/weights(ip2j);
    end
    if DATA.KIN_E
    % same for electrons!!
    [JJ,PP] = meshgrid(Je,Pe);
    P2Je = PP + 2*JJ;
    % form an axis of velocity ordered moments
    p2je = unique(P2Je);
    % weights to average
    weights = zeros(size(p2je));
    % space time of moments amplitude wrt to p+2j degree
    Ne_ST = zeros(numel(p2je),numel(DATA.Ts5D));
    % fill the st diagramm by averaging every p+2j=n moments together
    for ip = 1:numel(Pe)
        for ij = 1:numel(Je)
            [~,ip2j] = min(abs(p2ji-(Pe(ip)+2*Je(ij))));
            Ne_ST(ip2j,:) = Ne_ST(ip2j,:) + transpose(squeeze(Nepj(ip,ij,:)));
            weights(ip2j) = weights(ip2j) + 1;
        end
    end
    % doing the average
    for ip2j = 1:numel(p2ji)
        Ne_ST(ip2j,:) = Ne_ST(ip2j,:)/weights(ip2j);
    end 
    end
    % plots
    % scaling
%     plt = @(x,ip2j) x./max(max(x));
    if OPTIONS.NORMALIZED
    plt = @(x,ip2j) x(ip2j,:)./max(x(ip2j,:));
    else
    plt = @(x,ip2j) x;
    end
    if DATA.KIN_E
    subplot(2,1,1)
    end
        imagesc(DATA.Ts5D,p2ji,plt(Ni_ST,1:numel(p2ji))); 
        set(gca,'YDir','normal')        
%         pclr = pcolor(XX,YY,plt(Ni_ST));
%         set(pclr, 'edgecolor','none'); hold on;
    xlabel('$t$'); ylabel('$p+2j$')
    title('$\langle\sum_k |N_i^{pj}|\rangle_{p+2j=const}$')
    if DATA.KIN_E
    subplot(2,1,2)
        imagesc(DATA.Ts5D,p2je,plt(Ne_ST,1:numel(p2ji))); 
        set(gca,'YDir','normal')
%         pclr = pcolor(XX,YY,plt(Ne_ST)); 
%         set(pclr, 'edgecolor','none'); hold on;
    xlabel('$t$'); ylabel('$p+2j$')
    title('$\langle\sum_k |N_e^{pj}|\rangle_{p+2j=const}$')
    suptitle(DATA.param_title);
    end

end

end

