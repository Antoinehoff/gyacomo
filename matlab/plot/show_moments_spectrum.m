function [ FIGURE ] = show_moments_spectrum( DATA, OPTIONS )
species_name = {'i','e'}; % hard coded because this list is not output yet

Pa = DATA.grids.Parray;
Ja = DATA.grids.Jarray;
Time_ = DATA.Ts3D;
FIGURE.fig = figure; FIGURE.FIGNAME = ['mom_spectrum_',DATA.params_string];
% set(gcf, 'Position',  [100 50 1000 400])
if OPTIONS.ST == 0
    OPTIONS.LOGSCALE = 0;
end
if OPTIONS.LOGSCALE
    logname = 'log';
    compress = @(x,ia) log(sum(abs((x(:,:,:,:,:))),4));
else
    logname = '';
    compress = @(x,ia) sum(abs((x(:,:,:,:,:))),4);
end
for ia = 1:DATA.inputs.Na
    Napjz = compress(DATA.Napjz,ia);
    Napjz = reshape(Napjz,DATA.grids.Np,DATA.grids.Nj,numel(DATA.Ts3D)); 
    subplot(double(DATA.inputs.Na),1,double(ia))
    plotname = [logname,'$\langle\sum_k |N_',species_name{ia},'^{pj}|\rangle_{p+2j=const}$'];
    %We order the moments (0,0) (1,0) (2,0) (0,1) (3,0) (1,1) (4,0) (2,1) (0,2) etc.
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

    if OPTIONS.FILTER
    % Experimental filtering
    nt_fourth  = ceil(numel(Time_)/4);
    nt_half    = ceil(numel(Time_)/2);
    Na_ST_avg = mean(Na_ST(:,nt_fourth:nt_half),2);
    Na_ST     = (Na_ST - Na_ST_avg)./Na_ST;
    OPTIONS.NORMALIZED = 0;
    %
    end
    % plots
    % scaling
    if OPTIONS.TAVG_2D
        nt_half    = ceil(numel(Time_)/2);        
        Napjz_tavg = mean(Napjz(:,:,nt_half:end),3);
        x_ = DATA.grids.Parray;
        y_ = DATA.grids.Jarray;
        if OPTIONS.TAVG_2D_CTR
            [YY,XX] = meshgrid(y_,x_);
            surf(XX,YY,Napjz_tavg);
        else
            imagesc(x_,y_,Napjz_tavg');
            set(gca,'YDir','normal')        
            xlabel('$p$'); ylabel('$j$');
        end
        clb = colorbar; colormap(jet);
        clb.Label.String = '$\log\langle E(p,j) \rangle_t$';
        clb.TickLabelInterpreter = 'latex';
        clb.Label.Interpreter = 'latex';
        clb.Label.FontSize= 18;
    else
        if OPTIONS.NORMALIZED
        plt = @(x,i) reshape(x(i,:)./max(x(i,:)),numel(p2j),numel(Time_));
        else
        plt = @(x,i) reshape(x(i,:),numel(p2j),numel(Time_));
        end
        Nlabelmax = 15;
        nskip = floor(numel(p2j)/Nlabelmax);
        if OPTIONS.ST
            imagesc(Time_,p2j,plt(Na_ST,1:numel(p2j))); 
            set(gca,'YDir','normal')        
            xlabel('$t$'); ylabel(yname);
            if ~isempty(ticks_labels)
                yticks(p2j(1:nskip:end));
                yticklabels(ticks_labels(1:nskip:end))
            end
                if OPTIONS.FILTER
                caxis([0 0.2]);
                title('Relative fluctuation from the average');
                end
                colorbar
        else
            colors_ = jet(numel(p2j));
            for i = 1:numel(p2j)
               semilogy(Time_,squeeze(Na_ST(i,:)),...
                   'DisplayName',ticks_labels{i},...
                   'color',colors_(i,:)); hold on;
            end
        title(plotname)
        end
    end
top_title(DATA.paramshort)

end

