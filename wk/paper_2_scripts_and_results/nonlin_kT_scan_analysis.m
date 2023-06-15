kN=2.22;
figure
ERRBAR = 0; LOGSCALE = 1;
% nustr = '1e-3'; mrkstyl='o';
% nustr = '1e-2'; mrkstyl='^';
nustr = '5e-2'; mrkstyl='s';
GENE = 1;
if ~GENE
    rootdir = ['/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_',nustr];
    % mrkstyl='v';
else
    rootdir = '/misc/gene_results/kT_scan_nu0'; 
    mrkstyl='d';
end
Mmax = 6;
msz = 10; lwt = 2.0;
xname = '$\kappa_T (\kappa_N=2.22)$';
scanvarname = 'kT';
scanvalues = [9.0 8.0 6.96,6.5:-0.5:4.0];
tplotvalues= [9.0 8.0 6.96 6.5 6 5];
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
fclose(fid);
system('command rm list.txt');
[~,ids] = sort(Ps);
directories = directories(ids); Ps = Ps(ids);

if GENE
    % clrs_ = lines(numel(directories));
    clrs_ = cool(numel(directories));
else
    clrs_ = cool(numel(directories));
    % clrs_ = lines(numel(directories));
end
M   = numel(directories);
chi_kT_PJ = zeros(numel(scanvalues),M);
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
    itp = 1;
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
        try
            Trange  = data.Ts0D(end)*[0.5 1.0];
        catch % if data does not exist put 0 everywhere
            data.Ts0D = 0;
            data.HFLUX_X = 0;
            Trange = 0;
            data.inputs.PMAX = Ps(j);
            data.inputs.JMAX = Js(j);
            data.inputs.K_T  = kTs(i);
            data.inputs.K_N  = kN;
        end
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
        x(i) = kTs(i);
        if itp <= numel(tplotvalues)
        if scanvalues(i) == tplotvalues(itp)
            subplot(numel(tplotvalues),2,2*itp-1)
            hold on;
            Qx      = data.HFLUX_X;
            T       = data.Ts0D;
            plot(T,Qx,'DisplayName',...
                ['$Q_{avg}=',sprintf('%2.2f',Qx_avg(i)),'\pm',sprintf('%2.2f',Qx_std(i)),'$'],...
                'Color',clr_); hold on
            itp = itp + 1;
        end
        end
    end
    % plot;
    subplot(222)
    hold on;
    if ERRBAR
    errorbar(x,Chi_avg,Chi_std,'DisplayName',...
        ['(',num2str(data.inputs.PMAX),',',num2str(data.inputs.JMAX),')'],...
        'color',clr_,'Marker',mrkstyl,'MarkerFaceColor',clr_,...
        'MarkerSize',7,'LineWidth',2); 
    else
        plot(x,Chi_avg,'DisplayName',...
        ['(',num2str(data.inputs.PMAX),',',num2str(data.inputs.JMAX),')'],...
        'color',clr_,'Marker',mrkstyl,'MarkerFaceColor',clr_,...
        'MarkerSize',7,'LineWidth',2); 
    end
    hold on;
    chi_kT_PJ(1:N,j) = Chi_avg;
end
% Formatting and add GENE ref
for i = 1:numel(tplotvalues)
    subplot(numel(tplotvalues),2,2*i-1)
    ylabel('$Q_x$');
    yl = ylim; xl = xlim;
    title(['$R/L_T=',num2str(tplotvalues(i)),'$'],'Position',[xl(2)/4 yl(2)]);
    if LOGSCALE 
        set(gca,'YScale','log')
    else
        set(gca,'YScale','linear');
    end
    if i<numel(tplotvalues)
        xticklabels([]);
    else
        xlabel('$t c_s/R$');
    end
    grid off
    xlim([0 1000]);
end
subplot(222)
hold on;
Dim2000 = load('/home/ahoffman/gyacomo/wk/benchmark_and_scan_scripts/Dimits_2000_fig3_full_no_GF.txt');
plot(Dim2000(:,1),Dim2000(:,2),'ok','DisplayName','Dimits 2000',...
        'MarkerFaceColor','k','MarkerSize',7,'LineWidth',2);
% xline(4.0,'DisplayName','Dimits $\kappa_T^{crit}$','color',[0 0 0])
Dim2000 = load('/home/ahoffman/gyacomo/wk/benchmark_and_scan_scripts/Dimits_2000_fig3_dashed_low.txt');
plot(Dim2000(:,1),Dim2000(:,2),'--k','DisplayName','Dimits 2000 fit');
if LOGSCALE
    set(gca,'YScale','log')
end

ylabel('$\chi L_N/\rho_s^2 c_s$');
xlabel('$R/L_T$');
title(['$\nu_{DGDK}=$',nustr])
legend('show');
legend('Location','northwest')
xlim([3.5 10]);
grid on
% ylim([0 3.5]);

if 0
%%
% subplot(224)
figure
[PJ,KT] = meshgrid(Ps,scanvalues);
contourf(KT,PJ,(chi_kT_PJ),10)
% pclr=pcolor(KT,PJ,chi_kT_PJ);
xlabel('$R/L_T$')
ylabel('P (J=P/2)');
clb=colorbar; 
colormap(bluewhitered)
clb.Label.String = '$\chi$';
clb.Label.Interpreter = 'latex';
clb.Label.FontSize= 18;
end