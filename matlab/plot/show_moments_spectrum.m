function [ FIGURE ] = show_moments_spectrum( DATA, OPTIONS )

Pi = DATA.Pi;
Ji = DATA.Ji;
if ~(isempty(DATA.Nipjz))
    Time_ = DATA.Ts3D;
    Nipj = sum(abs(DATA.Nipjz),3);
else
    Time_ = DATA.Ts5D;
%     Nipj = sum(sum(sum(abs(DATA.Nipj).^2,3),4),5);
    Nipj = sum(sum(sum(conj(DATA.Nipj).*DATA.Nipj,3),4),5);
end
Nipj = squeeze(Nipj);

if DATA.KIN_E
Pe = DATA.Pe;
Je = DATA.Je;
    if ~(isempty(DATA.Nepjz))
        Nepj = sum(abs(DATA.Nepjz),3);
    else
        Nepj = sum(sum(sum(abs(DATA.Nepj),3),4),5);
    end
    Nepj = squeeze(Nepj);
end

FIGURE.fig = figure; FIGURE.FIGNAME = ['mom_spectrum_',DATA.PARAMS];
set(gcf, 'Position',  [100 50 1000 400])

if ~OPTIONS.ST
    t0 = OPTIONS.TIME(1);
    t1 = OPTIONS.TIME(end);
    [~,it0] = min(abs(t0-Time_));
    [~,it1] = min(abs(t1-Time_));
    Nipj = mean(Nipj(:,:,it0:it1),3);
    if DATA.KIN_E
    Nepj = mean(Nepj(:,:,it0:it1),3);
    end
    if numel(OPTIONS.TIME) == 1
        TITLE=['$t=',num2str(OPTIONS.TIME),'$'];
    else
        TITLE=['$t\in[',num2str(t0),',',num2str(t1),']$'];
    end 
    Nipj = squeeze(Nipj);

    ymini = min(Nipj); ymaxi = max(Nipj);

    if DATA.KIN_E
    Nepj = squeeze(Nepj);
    ymine = min(Nepj); ymaxe = max(Nepj);
    ymax  = max([ymaxi ymaxe]);
    ymin  = min([ymini ymine]);
    end
    if DATA.KIN_E
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
    
    if DATA.KIN_E
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
    if OPTIONS.P2J
        plotname = '$\langle\sum_k |N_i^{pj}|\rangle_{p+2j=const}$';
        [JJ,PP] = meshgrid(Ji,Pi);
        P2Ji = PP + 2*JJ;
        % form an axis of velocity ordered moments
        y_ = unique(P2Ji); yname = '$p+2j$';
        % weights to average
        weights = zeros(size(y_));
        % space time of moments amplitude wrt to p+2j degree
        Ni_ST = zeros(numel(y_),numel(Time_));
        % fill the st diagramm by averaging every p+2j=n moments together
        for ip = 1:numel(Pi)
            for ij = 1:numel(Ji)
                [~,ip2j] = min(abs(y_-(Pi(ip)+2*Ji(ij))));
                Ni_ST(ip2j,:) = Ni_ST(ip2j,:) + transpose(squeeze(Nipj(ip,ij,:)));
                weights(ip2j) = weights(ip2j) + 1;
            end
        end
        % doing the average
        for ip2j = 1:numel(y_)
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
        Ne_ST = zeros(numel(p2je),numel(Time_));
        % fill the st diagramm by averaging every p+2j=n moments together
        for ip = 1:numel(Pe)
            for ij = 1:numel(Je)
                [~,ip2j] = min(abs(y_-(Pe(ip)+2*Je(ij))));
                Ne_ST(ip2j,:) = Ne_ST(ip2j,:) + transpose(squeeze(Nepj(ip,ij,:)));
                weights(ip2j) = weights(ip2j) + 1;
            end
        end
        % doing the average
        for ip2j = 1:numel(y_)
            Ne_ST(ip2j,:) = Ne_ST(ip2j,:)/weights(ip2j);
        end 
        end
        ticks_labels = {};
    else % We just order our moments w.r.t. to the convention ..
         % (0,0) (1,0) (2,0) (0,1) (3,0) (1,1) (4,0) (2,1) (0,2) etc.
        plotname = '$\langle\sum_k |N_i^{pj}|^2\rangle_z$';
        Nmoments = numel(Nipj(:,:,1)); % total number of moments
        HL_deg   = zeros(Nmoments,2);   % list of degrees, first column Hermite, second Laguerre
        im  = 2;
        deg = 1; % the zero degree is always here first place
        ticks_labels = cell(10,1);
        ticks_labels{1} = '(0,0)';
        while(im<=Nmoments)
        FOUND = 1;
            while(FOUND) % As we find a pair matching the degree we retry
            FOUND = 0;
            for ij = 1:DATA.Jmaxi
            for ip = 1:DATA.Pmaxi
                if((ip-1)+2*(ij-1) == deg)
                    % Check if the pair is already added
                    check_ = HL_deg == [DATA.Pi(ip) DATA.Pi(ij)];
                    check_ = sum(check_(:,1) .* check_(:,2));
                    if ~check_ % if not add it
                        HL_deg(im,1) = DATA.Pi(ip);
                        HL_deg(im,2) = DATA.Pi(ij);
                        ticks_labels{im} = ['(',num2str(DATA.Pi(ip)),',',num2str(DATA.Ji(ij)),')'];
                        im = im + 1; FOUND = 1;
                    end
                end
                end
            end
            end
            % No pair found anymore, increase the degree
            deg = deg + 1;
        end
        
        % form an axis of velocity ordered moments
        y_ = 1:Nmoments; yname = '$(P,J)$';
        % space time of moments amplitude wrt to p+2j degree
        Ni_ST = zeros(Nmoments,numel(Time_));
        for i_ = 1:Nmoments
           Ni_ST(i_,:) = Nipj(HL_deg(i_,1)+1,HL_deg(i_,2)+1,:); 
        end
        if DATA.KIN_E
        % space time of moments amplitude wrt to p+2j degree
        Ne_ST = zeros(Nmoments,numel(Time_));
        for i_ = 1:Nmoments
           Ne_ST = Nepj(HL_deg(i_,1)+1,HL_deg(i_,2)+1,:); 
        end
        end    
         
    end
    % plots
    % scaling
    if OPTIONS.NORMALIZED
    plt = @(x,i) x(i,:)./max(x(i,:));
    else
    plt = @(x,i) x;
    end
    if DATA.KIN_E
    subplot(2,1,1)
    end
    
    imagesc(Time_,y_,plt(Ni_ST,1:numel(y_))); 
    set(gca,'YDir','normal')        
    xlabel('$t$'); ylabel(yname);
    if ~isempty(ticks_labels)
        yticks(y_);
        yticklabels(ticks_labels)
    end
    title(plotname)
    
    if DATA.KIN_E
    subplot(2,1,2)
        imagesc(Time_,p2je,plt(Ne_ST,1:numel(y_))); 
        set(gca,'YDir','normal')
        xlabel('$t$'); ylabel(yname)
        title(plotname)
        suptitle(DATA.param_title);
    end

end

end

