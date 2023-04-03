function [ FIGURE ] = show_moments_spectrum( DATA, OPTIONS )
species_name = {'i','e'}; % hard coded because this list is not output yet

Pa = DATA.grids.Parray;
Ja = DATA.grids.Jarray;
Time_ = DATA.Ts3D;
FIGURE.fig = figure; FIGURE.FIGNAME = ['mom_spectrum_',DATA.params_string];
set(gcf, 'Position',  [100 50 1000 400])
if OPTIONS.LOGSCALE
%     compress = @(x,ia) log(sum(abs(squeeze(x(ia,:,:,:))),3));
    compress = @(x,ia) log(sum(abs(squeeze(x(:,:,:,:))),3));
else
%     compress = @(x,ia) sum(abs(squeeze(x(ia,:,:,:))),3);
    compress = @(x,ia) sum(abs(squeeze(x(:,:,:,:))),3);
end
for ia = 1:DATA.inputs.Na
%     Napjz  = sum(abs(squeeze(DATA.Napjz(ia,:,:,:,:))),3);
    Napjz  =compress(DATA.Napjz);
    subplot(double(DATA.inputs.Na),1,double(ia))
    if OPTIONS.P2J
        plotname = ['$\langle\sum_k |N_',species_name{ia},'^{pj}|\rangle_{p+2j=const}$'];
        [JJ,PP] = meshgrid(Ja,Pa);
        P2J = PP + 2*JJ;
        % form an axis of velocity ordered moments
        p2j = unique(P2J); yname = '$p+2j$';
        % weights to average
        weights = zeros(size(p2j));
        % space time of moments amplitude wrt to p+2j degree
        Na_ST = zeros(numel(p2j),numel(Time_));
        % fill the st diagramm by averaging every p+2j=n moments together
        for ip = 1:numel(Pa)
            for ij = 1:numel(Ja)
                [~,ip2j] = min(abs(p2j-(Pa(ip)+2*Ja(ij))));
                Na_ST(ip2j,:) = Na_ST(ip2j,:) + transpose(squeeze(Napjz(ip,ij,:)));
                weights(ip2j) = weights(ip2j) + 1;
            end
        end
        % doing the average
        for ip2j = 1:numel(p2j)
            Na_ST(ip2j,:) = Na_ST(ip2j,:)/weights(ip2j);
        end
        ticks_labels = {};
    else % We just order our moments w.r.t. to the convention ..
         % (0,0) (1,0) (2,0) (0,1) (3,0) (1,1) (4,0) (2,1) (0,2) etc.
        plotname = '$\langle\sum_k |N_i^{pj}|^2\rangle_z$';
        Nmoments = numel(Napjz(:,:,1)); % total number of moments
        HL_deg   = zeros(Nmoments,2);   % list of degrees, first column Hermite, second Laguerre
        im  = 2;
        deg = 1; % the zero degree is always here first place
        ticks_labels = cell(10,1);
        ticks_labels{1} = '(0,0)';
        while(im<=Nmoments)
        FOUND = 1;
            while(FOUND) % As we find a pair matching the degree we retry
            FOUND = 0;
            for ij = 1:DATA.grids.Nj
                for ip = 1:DATA.grids.Np
                    if((ip-1)+2*(ij-1) == deg)
                        % Check if the pair is already added
                        check_ = HL_deg == [DATA.grids.Parray(ip) DATA.grids.Jarray(ij)];
                        check_ = sum(check_(:,1) .* check_(:,2));
                        if ~check_ % if not add it
                            HL_deg(im,1) = DATA.grids.Parray(ip);
                            HL_deg(im,2) = DATA.grids.Jarray(ij);
                            ticks_labels{im} = ['(',num2str(DATA.grids.Parray(ip)),',',num2str(DATA.grids.Jarray(ij)),')'];
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
        p2j = 1:Nmoments; yname = '$(P,J)$';
        % space time of moments amplitude wrt to p+2j degree
        Na_ST = zeros(Nmoments,numel(Time_));
        for i_ = 1:Nmoments
           Na_ST(i_,:) = Napjz((HL_deg(i_,1)/DATA.grids.deltap)+1,HL_deg(i_,2)+1,:); 
        end
    end    
    % plots
    % scaling
    if OPTIONS.NORMALIZED
    plt = @(x,i) x(i,:)./max(x(i,:));
    else
    plt = @(x,i) x;
    end

    imagesc(Time_,p2j,plt(Na_ST,1:numel(p2j))); 
    set(gca,'YDir','normal')        
    xlabel('$t$'); ylabel(yname);
    if ~isempty(ticks_labels)
        yticks(p2j);
        yticklabels(ticks_labels)
    end
    title(plotname)

end
suptitle(DATA.paramshort)

end

