function [ FIGURE ] = show_napjz( DATA, OPTIONS )
fname = DATA.outfilenames{OPTIONS.JOBNUM+1};
switch OPTIONS.specie
    case 'i'
        Pa = DATA.Pi;
        Ja = DATA.Ji;
        Napjz = load_pjz_data(fname,'Nipjz','i');
        name = 'N_i^{pj}';
        FIGURE.FIGNAME = ['Nipj_spectrum_',DATA.PARAMS];
    case 'e'
        Pa = DATA.Pe;
        Ja = DATA.Je;
        Napjz = load_pjz_data(fname,'Nipjz','i');
        name = 'N_e^{pj}';
        FIGURE.FIGNAME = ['Nepj_spectrum_',DATA.PARAMS];
end

switch OPTIONS.compz
    case 'avg'
        Napjz = squeeze(mean(Napjz,3));
    otherwise
        Napjz = squeeze(Napjz(:,:,OPTIONS.compz));
end

FIGURE.fig = figure; 
set(gcf, 'Position',  [100 50 1000 400])

switch OPTIONS.PLOT_TYPE
    case 'space-time'
    [JJ,PP] = meshgrid(Ja,Pa);
    P2Ja = PP + 2*JJ;
    % form an axis of velocity ordered moments
    p2ja = unique(P2Ja);
    % weights to average
    weights = zeros(size(p2ja));
    % space time of moments amplitude wrt to p+2j degree
    Na_ST = zeros(numel(p2ja),numel(DATA.Ts5D));
    % fill the st diagramm by averaging every p+2j=n moments together
    for ip = 1:numel(Pa)
        for ij = 1:numel(Ja)
            [~,ip2j] = min(abs(p2ja-(Pa(ip)+2*Ja(ij))));
            Na_ST(ip2j,:) = Na_ST(ip2j,:) + transpose(squeeze(Napjz(ip,ij,:)));
            weights(ip2j) = weights(ip2j) + 1;
        end
    end
    % doing the average
    for ip2j = 1:numel(p2ja)
        Na_ST(ip2j,:) = Na_ST(ip2j,:)/weights(ip2j);
    end
    % plots
    if OPTIONS.NORMALIZED
    plt = @(x,ip2j) x(ip2j,:)./max(x(ip2j,:));
    else
    plt = @(x,ip2j) x;
    end
    imagesc(DATA.Ts5D,p2ja,plt(Na_ST,1:numel(p2ja))); 
    set(gca,'YDir','normal')        
    xlabel('$t$'); ylabel('$p+2j$')
    title('$\langle\sum_k |',name,'|\rangle_{p+2j=const}$')
    if DATA.K_E
    subplot(2,1,2)
        imagesc(DATA.Ts5D,p2je,plt(Ne_ST,1:numel(p2ja))); 
        set(gca,'YDir','normal')
    xlabel('$t$'); ylabel('$p+2j$')
    title('$\langle\sum_k |N_e^{pj}|\rangle_{p+2j=const}$')
    suptitle(DATA.param_title);
    end
    
    case 'Tavg-1D'
    t0 = OPTIONS.TIME(1);
    t1 = OPTIONS.TIME(end);
    [~,it0] = min(abs(t0-DATA.Ts5D));
    [~,it1] = min(abs(t1-DATA.Ts5D));
    Napjz = mean(Napjz(:,:,it0:it1),3);
    if DATA.K_E
    Napjz = mean(Napjz(:,:,it0:it1),3);
    end
    if numel(OPTIONS.TIME) == 1
        TITLE=['$t=',num2str(OPTIONS.TIME),'$'];
    else
        TITLE=['$t\in[',num2str(t0),',',num2str(t1),']$'];
    end 
    Napjz = squeeze(Napjz);

    ymini = min(Napjz); ymaxi = max(Napjz);

    if DATA.K_E
    Napjz = squeeze(Napjz);
    ymine = min(Napjz); ymaxe = max(Napjz);
    ymax  = max([ymaxi ymaxe]);
    ymin  = min([ymini ymine]);
    end
    if DATA.K_E
    subplot(121)
    end
    if ~OPTIONS.P2J
        for ij = 1:numel(Ja)
            name = ['$j=',num2str(Ja(ij)),'$'];
            semilogy(Pa,Napjz(:,ij),'o-','DisplayName',name); hold on;
        end
        xlabel('$p$'); 
    else
        for ij = 1:numel(Ja)
            name = ['$j=',num2str(Ja(ij)),'$'];
            semilogy(Pa+2*Ja(ij),Napjz(:,ij),'o-','DisplayName',name); hold on;
        end
        xlabel('$p+2j$')
    end
    ylabel(['$\sum_{kx,ky}|N_i^{pj}|$']);
    legend('show');
    title([TITLE,' He-La ion spectrum']);
    
    if DATA.K_E
    subplot(122)
    if ~OPTIONS.P2J
        for ij = 1:numel(Ja)
            name = ['$j=',num2str(Ja(ij)),'$'];
            semilogy(Pa,Napjz(:,ij),'o-','DisplayName',name); hold on;
        end
        xlabel('p'); 
    else
    for ij = 1:numel(Ja)
        name = ['$j=',num2str(Ja(ij)),'$'];
        semilogy(Pa+2*Ja(ij),Napjz(:,ij),'o-','DisplayName',name); hold on;
    end
    xlabel('$p+2j$')
    end
    ylabel(['$\sum_{kx,ky}|N_e^{pj}|$']);
    legend('show');
    title([TITLE,' He-La elec. spectrum']);
    end

end

end

