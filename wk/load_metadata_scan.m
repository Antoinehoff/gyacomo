%% Metadata locations
% Scans over KT and PJ, keeping ky, CO constant
% datafname = 'p2_CBC_convergence_KT_PJ/12x24_ky_0.3_kT_3_7__P_2_30_DGdkaa_0.05.mat';
% datafname = 'p2_CBC_convergence_KT_PJ/12x24_ky_0.3_kT_3_7__P_2_30_DGdkaa_0.01.mat';
% datafname = 'p2_CBC_convergence_KT_PJ/12x24_ky_0.3_kT_3_7__P_2_30_DGdkaa_0.mat';
% Scans over NU and PJ, keeping ky and KY constant
% datafname = 'p2_CBC_convergence_coll_PJ/12x24_ky_0.3_nu_0_0.1__P_2_30_KT_6.96.mat';
datafname = 'p2_CBC_convergence_coll_PJ/12x24_ky_0.3_nu_0_0.1__P_2_30_KT_5.3.mat';

%% Load data
fname = ['../results/',datafname];
d = load(fname);
if 1
%% Pcolor of the peak
figure;
[XX_,YY_] = meshgrid(d.s1,d.s2);
pclr=imagesc_custom(XX_,YY_,d.data'.*(d.data>0)');
title(d.title);
xlabel(d.s1name); ylabel(d.s2name);
colormap(bluewhitered)
clb=colorbar; 
clb.Label.String = '$\gamma c_s/R$';
clb.Label.Interpreter = 'latex';
clb.Label.FontSize= 18;
end
if 0
%% Scan along first dimension
figure
colors_  = jet(numel(d.s2));
for i = 1:numel(d.s2)
%     plot(d.s1,d.data(:,i),'s-',...
%         'LineWidth',2.0,...
%         'DisplayName',[d.s2name,'=',num2str(d.s2(i))],...
%         'color',colors_(i,:)); 
    errorbar(d.s1,d.data(:,i),d.err(:,i),'s-',...
    'LineWidth',2.0,...
    'DisplayName',[d.s2name,'=',num2str(d.s2(i))],...
    'color',colors_(i,:)); 
    hold on;
end
xlabel(d.s1name); ylabel(d.dname);title(d.title);
xlim([d.s1(1) d.s1(end)]);
colormap(colors_);
clb = colorbar;
caxis([d.s2(1),d.s2(end)]);
clb.YTick=d.s2;
clb.Label.String = d.s2name;
clb.Label.Interpreter = 'latex';
clb.Label.FontSize= 18;
end
if 0
%% Scan along second dimension
figure
colors_ = jet(numel(d.s1));
for i = 1:numel(d.s1)
%     plot(d.s2,d.data(i,:),'s-',...
%         'LineWidth',2.0,...
%         'DisplayName',[d.s1name,'=',num2str(d.s1(i))],...
%         'color',colors_(i,:)); 
    errorbar(d.s2,d.data(i,:),d.err(i,:),'s-',...
        'LineWidth',2.0,...
        'DisplayName',[d.s1name,'=',num2str(d.s1(i))],...
        'color',colors_(i,:)); 
    hold on;
end
xlabel(d.s2name); ylabel(d.dname);title(d.title);
xlim([d.s2(1) d.s2(end)]);
colormap(jet);
clb = colorbar;
caxis([d.s1(1),d.s1(end)]);
clb.YTick=d.s1;
clb.Label.String = d.s1name;
clb.Label.Interpreter = 'latex';
clb.Label.FontSize= 18;
end

if 0
%% Convergence analysis
figure
colors_ = jet(numel(d.s1));
for i = 1:numel(d.s1)
    target_ = d.data(i,end);
    eps_    = abs(target_ - d.data(i,1:end-1));
    semilogy(d.s2(1:end-1),eps_,'s',...
        'LineWidth',2.0,...
        'DisplayName',[d.s1name,'=',num2str(d.s1(i))],...
        'color',colors_(i,:));
    hold on;
end
xlabel(d.s2name); ylabel('$\epsilon$');title(d.title);
xlim([d.s2(1) d.s2(end)]);
colormap(jet);
clb = colorbar;
caxis([d.s1(1),d.s1(end)]);
clb.YTick=d.s1;
clb.Label.String = d.s1name;
clb.Label.Interpreter = 'latex';
clb.Label.FontSize= 18;
grid on;
end