kN=2.22;
figure
ERRBAR = 0; LOGSCALE = 0;
nustr = '1e-3';
rootdir = ['/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_',nustr];GENE = 0;
% rootdir = '/misc/gene_results/kT_scan_nu0'; GENE = 1;
msz = 10; lwt = 2.0;
mrkstyl='v';
xname = '$\kappa_T (\kappa_N=2.22)$';
scanvarname = 'kT';
scanvalues = [6.96,6.5:-0.5:4.0];
% Get all subdirectories
system(['ls -d ',rootdir,'/*/ > list.txt']);
fid = fopen('list.txt');
tline = fgetl(fid); i_ = 1; Ps=[]; Js =[]; directories={};
while ischar(tline)
    directories{i_} =  tline;
    resstr{i_} = tline(numel(rootdir)+2:end-1);
    tmp = sscanf(resstr{i_},'%dx%dx%dx%dx%d');
    Ps  = [Ps tmp(1)];
    Js  = [Js tmp(2)];
    tline = fgetl(fid);
    i_ = i_+1;
end
[~,ids] = sort(Ps);
fclose(fid);
system('command rm list.txt');5

directories = directories(ids); Ps = Ps(ids);
if GENE
    clrs_ = jet(numel(directories));
else
    clrs_ = cool(numel(directories));
end
M = numel(directories);

% chi_kT_PJ = zeros(numel(scanvalues),M);
for j = 1:M
     % Get all subdirectories
    system(['ls -d ',directories{j},'*/ > list.txt']);
    fid = fopen('list.txt');
    tline = fgetl(fid); i_ = 1; kTs=[]; subdirectories={};
    while ischar(tline)
        subdirectories{i_} =  tline;
        str = tline(numel(directories{j})+1:end-1);
        tmp = sscanf(str,'kT_%f');
        kTs  = [kTs tmp(1)];
        tline = fgetl(fid);
        i_ = i_+1;
    end
    fclose(fid);
    system('command rm list.txt');
    [~,ids] = sort(kTs,'descend');
    subdirectories = subdirectories(ids); kTs = kTs(ids);

    naming = @(s) sprintf('%1.1f',s); clr_ = clrs_(j,:);
    N   = numel(subdirectories);
    x = 1:N;
    Qx_avg  = 1:N;
    Qx_std  = 1:N;
    Chi_avg = 1:N;
    Chi_std = 1:N;
    data = {};
    for i = 1:N
        subdir = subdirectories{i};
        if ~GENE
            data    = compile_results_low_mem(data,subdir,00,10);
        else
            namelist      = read_namelist([subdir,'parameters']);
            data.inputs.PMAX = namelist.box.nv0;
            data.inputs.JMAX = namelist.box.nw0;
            data.inputs.K_T  = kTs(i);
            data.inputs.K_N  = namelist.species.omn;
            nrgfile          = 'nrg.dat.h5';
            data.Ts0D        = h5read([subdir,nrgfile],'/nrgions/time');
            data.HFLUX_X     = h5read([subdir,nrgfile],'/nrgions/Q_es');
        end
        Trange  = data.Ts0D(end)*[0.3 1.0];
        %
        [~,it0] = min(abs(Trange(1)  -data.Ts0D)); 
        [~,it1] = min(abs(Trange(end)-data.Ts0D)); 
        %
        if 0
            Qxa_    = 0*(1:Nseg);
            for n = 1:Nseg
               ntseg = floor((it1-it0)/n);
               for m = 1:n 
                Qxa_(n) = Qxa_(n) + mean(Qx((1:ntseg)+(m-1)*ntseg));
               end
               Qxa_(n) = Qxa_(n)/n;
            end
            Qx_avg(i) = mean(Qxa_);
            Qx_std(i) = std(Qxa_);
        else
            Qx_avg(i) = mean(data.HFLUX_X(it0:it1));
            Qx_std(i) =  std(data.HFLUX_X(it0:it1));
        end
        Chi_avg(i) = Qx_avg(i)./data.inputs.K_T/data.inputs.K_N;
        Chi_std(i) = Qx_std(i)./data.inputs.K_T/data.inputs.K_N;
        x(i) = data.inputs.K_T;
        subplot(N,2,2*i-1)
        hold on;
        Qx      = data.HFLUX_X;
        T       = data.Ts0D;
        plot(T,Qx,'DisplayName',...
            ['$Q_{avg}=',sprintf('%2.2f',Qx_avg(i)),'\pm',sprintf('%2.2f',Qx_std(i)),'$'],...
            'Color',clr_); hold on
    end
    % plot;
    subplot(222)
    hold on;
    if ERRBAR
    errorbar(x,Chi_avg,Chi_std,'DisplayName',...
        ['(',num2str(data.inputs.PMAX),',',num2str(data.inputs.JMAX),')'],...
        'color',clr_,'Marker',mrkstyl); 
    else
        plot(x,Chi_avg,'DisplayName',...
        ['(',num2str(data.inputs.PMAX),',',num2str(data.inputs.JMAX),')'],...
        'color',clr_,'Marker',mrkstyl); 
    end
    hold on;
    % chi_kT_PJ(:,j) = Chi_avg;
end
% Formatting
for i = 1:N
    subplot(N,2,2*i-1)
    ylabel('$Q_x$');
    yl = ylim; xl = xlim;
    title(['$\kappa_T=',num2str(x(i)),'$'],'Position',[xl(2)/2 yl(2)]);
    if LOGSCALE 
        set(gca,'YScale','log')
    else
        set(gca,'YScale','linear');
    end
    if i<N
        xticklabels([]);
    else
        xlabel('$t c_s/R$');
    end
end
subplot(222)
hold on;
Dim2000 = load('/home/ahoffman/gyacomo/wk/benchmark_and_scan_scripts/Dimits_2000_fig3_full_no_GF.txt');
plot(Dim2000(:,1),Dim2000(:,2),'ok','DisplayName','Dimits 2000');
xline(4.0,'DisplayName','Dimits $\kappa_T^{crit}$','color',[0 0 0])
if LOGSCALE
    set(gca,'YScale','log')
end
	%-------------- GENE ---------------
kT_Qi_GENE = ...
    [...
     13. 2.7e+2 2.2e+1;...%128x64x16x24x12 kymin=0.02 (large box)
     11. 1.9e+2 1.7e+1;...%128x64x16x24x12 kymin=0.02 (large box)
     9.0 1.1e+2 4.2e+1;...%128x64x16x24x12 kymin=0.05
     7.0 3.5e+1 4.6e+0;...%128x64x16x24x12 kymin=0.05
     5.3 9.7e+0 6.8e+0;...%128x64x16x24x12 kymin=0.05
     4.5 2.3e-1 5.0e-2;...%128x64x16x24x12 kymin=0.05
    ];
y_ = kT_Qi_GENE  (:,2)./kT_Qi_GENE  (:,1)./kN;
e_ = kT_Qi_GENE  (:,3)./kT_Qi_GENE  (:,1)./kN;
plot(kT_Qi_GENE (:,1),y_,...
    '+-.k','DisplayName','GENE 24x12',...
    'MarkerSize',msz,'LineWidth',lwt); hold on
kT_Qi_GENE = ...
    [...
     6.5 1.9e+0;...%128x64x16x32x16 kymin=0.05
     6.0 1.0e-1;...%128x64x16x32x16 kymin=0.05
     5.5 3.9e-2;...%128x64x16x32x16 kymin=0.05
     5.3 1.4e-2;...%128x64x16x32x16 kymin=0.05
     4.5 1.2e-2;...%128x64x16x32x16 kymin=0.05
     4.0 2.9e-6;...%128x64x16x32x16 kymin=0.05
    ];
y_ = kT_Qi_GENE  (:,2);
plot(kT_Qi_GENE (:,1),y_,...
    '*-.k','DisplayName','GENE 32x16',...
    'MarkerSize',msz,'LineWidth',lwt); hold on
ylabel('$\chi$');
xlabel(xname);
title(['$\nu_{DGDK}=$',nustr])
legend('show');
legend('Location','northwest')
xlim([3.5 8]);
ylim([0 3.5]);

%%
% subplot(224)
% figure
% clrs_ = cool(N);
% % [PJ,KT] = meshgrid(Ps,scanvalues);
% % surf(KT,PJ,chi_kT_PJ)
% for i=1:N
%     target = chi_kT_PJ(i,end);
%     % loglog(Ps,abs(chi_kT_PJ(i,:)-target)/target,'o--','color',clrs_(i,:),...
%         % 'DisplayName',['$\kappa_T=$',num2str(x(i))]); hold on
%     semilogy(Ps,abs(chi_kT_PJ(i,:)-target)/target,'o--','color',clrs_(i,:),...
%         'DisplayName',['$\kappa_T=$',num2str(x(i))]); hold on
% end
% % ylabel('error in \%')
% ylabel('$\chi$')
% xlabel('P (J=P/2)');
