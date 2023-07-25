kN=2.22;
figure
ERRBAR = 1; LOGSCALE = 0; AU = 0;
resstr={};
msz = 10; lwt = 2.0;
% CO = 'DGGK'; mrkstyl='d';
% CO = 'SGGK'; mrkstyl='s'; 
CO = 'LDGK'; mrkstyl='o';
% GRAD = 'kT_7.0'; kT = 6.96;
GRAD = 'kT_5.3'; kT = 5.3;
% GRAD = 'kT_4.5'; kT = 4.5;
xname = ['$\nu_{',CO,'}$ '];
titlename = [CO,', ',GRAD];
scanvarname = 'nu';
rootdir = ['/misc/gyacomo23_outputs/paper_2_GYAC23/collision_study/nu',CO,'_scan_',GRAD]; 
scanval = {'0.005' '0.01' '0.02' '0.05' '0.1' '0.2' '0.5'};
naming = @(s) num2str(s);

 % Get all directories
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
system('command rm list.txt');

directories = directories(ids); Ps = Ps(ids);
% clrs_ = cool(numel(directories));
clrs_ = lines(numel(directories));

for j = 1:numel(directories)
     % Get all subdirectories
    system(['ls -d ',directories{j},'*/ > list.txt']);
    fid = fopen('list.txt');
    tline = fgetl(fid); i_ = 1; nus=[]; subdirectories={};
    while ischar(tline)
        subdirectories{i_} =  tline;
        str = tline(numel(directories{j})+1:end-1);
        tmp = sscanf(str,'nu_%f');
        nus  = [nus tmp(1)];
        tline = fgetl(fid);
        i_ = i_+1;
    end
    fclose(fid);
    system('command rm list.txt');
    [~,ids] = sort(nus,'descend');
    subdirectories = subdirectories(ids); nus = nus(ids);

    naming = @(s) sprintf('%1.1f',s); clr_ = clrs_(j,:);
    titlename = [CO,', ',GRAD,', ',resstr{j}];
    N   = numel(subdirectories);
    x = 1:N;
    Qx_avg  = 1:N;
    Qx_std  = 1:N;
    Chi_avg = 1:N;
    Chi_std = 1:N;
    data = {};
    for i = 1:N
        subdir = subdirectories{i};
        data    = compile_results_low_mem(data,subdir,00,20);
        try
            Trange  = data.Ts0D(end)*[0.5 0.75];
        catch % if data does not exist put 0 everywhere
            data.Ts0D = 0;
            data.HFLUX_X = 0;
            Trange = 0;
            data.inputs.PMAX = Ps(j);
            data.inputs.JMAX = Js(j);
            data.inputs.K_T  = kT;
            data.inputs.K_N  = kN;
            data.inputs.NU   = nus(i);
        end
            Trange  = data.Ts0D(end)*[0.25 0.75];
            % Trange  = [200 400];
        %
        [~,it0] = min(abs(Trange(1)  -data.Ts0D)); 
        [~,it1] = min(abs(Trange(end)-data.Ts0D)); 
        %
        Qx_avg(i) = mean(data.HFLUX_X(it0:it1));
        Qx_std(i) =  std(data.HFLUX_X(it0:it1));
        Chi_avg(i) = Qx_avg(i)./data.inputs.K_T/data.inputs.K_N;
        Chi_std(i) = Qx_std(i)./data.inputs.K_T/data.inputs.K_N;
        x(i) = data.inputs.NU;
        subplot(N,2,2*i-1)
        hold on;
        Qx      = data.HFLUX_X;
        if AU 
            Qx = Qx./max(Qx); 
        end
        T       = data.Ts0D;
        % Plot heatflux vs time
        plot(T,Qx,'DisplayName',[scanvarname,'=',num2str(x(i))],...
            'Color',clr_); hold on
        % plot([T(it0) T(end)],Qx_avg(i)*[1 1],'--k','DisplayName',...
        % ['$Q_{avg}=',sprintf('%2.2f',Qx_avg(i)),'\pm',sprintf('%2.2f',Qx_std(i)),'$']);
    end
    % plot chi vs nu
    subplot(122)
    hold on;
    if ERRBAR
    errorbar(x,Chi_avg,Chi_std,'DisplayName',...
        ['(',num2str(data.inputs.PMAX),',',num2str(data.inputs.JMAX),')'],...
        'color',clr_,'MarkerFaceColor','k','MarkerSize',7,'LineWidth',2); hold on;
    else
        plot(x,Chi_avg,'DisplayName',...
        ['(',num2str(data.inputs.PMAX),',',num2str(data.inputs.JMAX),')'],...
        'color',clr_,'Marker',mrkstyl); 
    end    
end

% Formatting
for i = 1:N
    subplot(N,2,2*i-1)
    ylabel('$Q_x$');
    yl = ylim; xl = xlim;
    yl(1) = 0; ylim(yl);
    title(['$\nu =',num2str(nus(i)),'$'],'Position',[xl(2)/2 yl(2)]);
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

% ------------- LIN kT=5.3 results
subplot(122)
hold on;
Lin1999 = load('/home/ahoffman/gyacomo/wk/benchmark_and_scan_scripts/Lin_1999_fig2.txt');
epsilon   = 0.18; q0 = 1.4;
nuLINvsGM = 0.5/epsilon^(3/2)*q0;
plot(Lin1999(:,1)/nuLINvsGM,Lin1999(:,2),'--ok','DisplayName','Lin1999');
yline(0.035348,'--k','DisplayName','$\nu\sim 0$')
set(gca,'XScale','log')
xlim([ 1e-3 1])
ylabel('$\chi$');
xlabel(xname);
title(titlename)
legend('show');
