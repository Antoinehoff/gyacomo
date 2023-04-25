clrs_ = lines(10);
kN=2.22;

% resdir = '5x3x128x64x24'; clr_ = clrs_(1,:);   CBC_res = [44.08 06.51]; kTthresh = 2.8;
% resdir = '7x4x128x64x24'; clr_ = clrs_(2,:);   CBC_res = [44.08 06.51]; kTthresh = 2.9;
% resdir = '9x2x128x64x24'; clr_ = clrs_(3,:); CBC_res = [37.60 04.65]; kTthresh = 3.3;
% resdir = '9x2x128x64x24_Lx200'; clr_ = clrs_(3,:); CBC_res = [37.60 04.65]; kTthresh = 3.3;
% resdir = '9x5x128x64x24'; clr_ = clrs_(3,:); CBC_res = [37.60 04.65]; kTthresh = 3.3;
% resdir = '9x5x128x64x24_Lx200'; clr_ = clrs_(3,:); CBC_res = [37.60 04.65]; kTthresh = 3.3;
% resdir = '13x2x128x64x24'; clr_ = clrs_(4,:); CBC_res = [36.84 07.22]; kTthresh = 4.4;
resdir = '13x5x128x64x24'; clr_ = clrs_(4,:); CBC_res = [37.60 04.65]; kTthresh = 3.9;

if 1
    scantype = 'kTscan';
    xname = '$\kappa_T (\kappa_N=2.22)$';
    scanvarname = 'kT';
    GRAD = '';
    prefix = ['/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_1e-3/',resdir,'/']; 
    scanval = [6.96 6.5 6.0 5.5 5.0 4.5 4.0];
    naming = @(s) sprintf('%1.1f',s);
    titlename = [GRAD,' $\nu=0$'];
else
    CO = 'DGGK'; mrkstyl='v'; clr_ = clrs_(2,:);
    % CO = 'SGGK'; mrkstyl='s'; clr_ = clrs_(1,:);
    % CO = 'LDGK'; mrkstyl='d'; clr_ = clrs_(5,:);
    % GRAD = 'CBC';
    GRAD = 'kT_5.3';
    % GRAD = 'kT_4.5';
    scantype = 'nuscan';
    xname = ['$\nu_{',CO,'}$ '];
    titlename = [CO,', ',resdir,', ',GRAD];
    scanvarname = 'nu';
    prefix = ['/misc/gyacomo23_outputs/paper_2_GYAC23/collision_study/nu',CO,'_scan_',GRAD,'/',resdir,'/']; 
    scanval = [0.005 0.01 0.02 0.05 0.1 0.2 0.5];
    naming = @(s) num2str(s);
end

N   = numel(scanval);
x = 1:N;
Qx_avg  = 1:N;
Qx_std  = 1:N;
Chi_avg = 1:N;
Chi_std = 1:N;
data = {};
figure

 
for i = 1:N
    datadir = [prefix,scanvarname,'_',naming(scanval(i)),'/'];
    Nseg = 5;

    data    = compile_results_low_mem(data,datadir,00,10);
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
    switch scantype
        case 'nuscan'
        x(i) = data.inputs.NU;
        case 'kTscan'
        x(i) = data.inputs.K_T;
    end
    subplot(N,2,2*i-1)
    Qx      = data.HFLUX_X;
    T       = data.Ts0D;
    plot(T,Qx,'DisplayName',[scanvarname,'=',num2str(x(i))],...
        'Color',clr_); hold on
    plot([T(it0) T(end)],Qx_avg(i)*[1 1],'--k','DisplayName',...
    ['$Q_{avg}=',sprintf('%2.2f',Qx_avg(i)),'\pm',sprintf('%2.2f',Qx_std(i)),'$']);
    ylabel(['$Q_x, (\kappa_T=',num2str(x(i)),')$']);
    % legend('show')
end
% Add colless CBC results
switch  scantype
    case 'kTscan'
        Chi_avg = [CBC_res(1)/2.22/6.96 Chi_avg];
        Chi_std = [CBC_res(2)/2.22/6.96 Chi_std];
        x = [6.96 x];
    case 'nuscan'
        Chi_avg = [CBC_res(1)/2.22/6.96 Chi_avg];
        Chi_std = [CBC_res(2)/2.22/6.96 Chi_std];
        x = [0 x];
end

% plot;
subplot(122)
errorbar(x,Chi_avg,Chi_std,'DisplayName',data.paramshort,'color',clr_,'Marker',mrkstyl); hold on;
switch  scantype
    case 'nuscan'
    Lin1999 = load('/home/ahoffman/gyacomo/wk/benchmark_and_scan_scripts/Lin_1999_fig2.txt');
    plot(Lin1999(:,1),Lin1999(:,2),'--ok','DisplayName','Lin1999');
    set(gca,'XScale','log')
    xlim([min(x)*0.8 max(x)*1.2])
    case 'kTscan'
    Dim2000 = load('/home/ahoffman/gyacomo/wk/benchmark_and_scan_scripts/Dimits_2000_fig3_full_no_GF.txt');
    plot(Dim2000(:,1),Dim2000(:,2),'ok','DisplayName','Dimits 2000');
    xline(kTthresh,'DisplayName','$\kappa_T^{crit}$','color',clr_)
    xline(4.0,'DisplayName','Dimits $\kappa_T^{crit}$','color',[0 0 0])
    xlim([kTthresh*0.8 max(x)*1.2])
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
	errorbar(kT_Qi_GENE (:,1),y_, e_,...
        '+-.k','DisplayName','GENE 128x64x16x24x12',...
	    'MarkerSize',msz,'LineWidth',lwt); hold on
end
% plot(ITG_threshold*[1 1],[0 20],'-.','DisplayName','$\kappa_T^{crit}$',...
%     'color',clr_);
ylabel('$\chi$');
xlabel(xname);
title(titlename)
% ylim(ylimits);
legend('show');
