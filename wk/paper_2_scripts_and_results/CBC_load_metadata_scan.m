% Metadata path
gyacomodir  = pwd;
gyacomodir = gyacomodir(1:end-2);
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add% EXECNAME = 'gyacomo_1.0';
%% Scans over KT and NU, keeping ky, CO constant (CBC_kT_nu_scan.m)
% datafname = 'p2_linear/8x24_ky_0.3_P_4_J_2_kT_3_6.96_nu_0_1_DGDK.mat';
% datafname = 'p2_linear/8x24_ky_0.3_P_8_J_4_kT_3_6.96_nu_0_1_DGDK.mat';
% datafname = 'p2_linear/8x24_ky_0.3_P_16_J_8_kT_3_6.96_nu_0_1_DGDK.mat';
%% Scans over KT and PJ, keeping ky, CO constant (CBC_kT_PJ_scan.m)
% datafname = 'p2_linear_new/8x24_ky_0.3_kT_3_6.96_P_2_30_DGDK_0.001.mat';
% datafname = 'p2_linear_check/8x24_ky_0.3_kT_3_9_P_2_30_DGDK_0.001.mat';
% datafname = 'p2_linear/8x24_ky_0.3_kT_3_6.96_P_2_40_DGDK_0.05.mat';
% datafname = 'p2_linear/8x24_ky_0.3_kT_3_6.96_P_2_30_DGDK_0.05.mat';
% datafname = 'p2_linear/8x24_ky_0.3_kT_3_6.96_P_2_40_DGDK_0.025.mat';
%% Scans over NU and PJ, keeping ky and KY constant (CBC_nu_PJ_scan.m)
% datafname = 'p2_linear/8x24_ky_0.3_nu_0.01_1_DGDK_P_2_40_KT_6.96.mat';
% % datafname = 'p2_linear/8x24_ky_0.3_nu_0.01_1_DGDK_P_2_40_KT_5.3.mat';
% datafname = 'p2_linear/8x24_ky_0.3_nu_0.01_1_SGGK_P_2_20_KT_6.96.mat';
 % datafname = 'p2_linear/8x24_ky_0.3_nu_0.01_1_hypcoll_P_2_30_KT_6.96.mat';
%  datafname = 'p2_linear/8x24_ky_0.3_nu_0.01_1_hypcoll_P_2_30_KT_5.3.mat';
 % datafname = 'p2_linear/8x24_ky_0.3_nu_0.01_1_dvpar4_P_2_30_KT_6.96.mat';
%  datafname = 'p2_linear/8x24_ky_0.3_nu_0.01_1_dvpar4_P_2_10_KT_6.96.mat';
%  datafname = 'p2_linear/8x24_ky_0.3_nu_0.01_1_dvpar4_P_6_10_KT_6.96.mat';
%  datafname = 'p2_linear/8x24_ky_0.3_nu_0.01_1_dvpar4_P_2_10_KT_5.3.mat';
%% Scans over P and J, keeping nu, ky and kT constant (CBC_P_J_scan.m)
% datafname = 'p2_linear/8x24_ky_0.3_J_2_15_P_2_30_kT_6.96_DGDK_0.05.mat';
% datafname = 'p2_linear/8x24_ky_0.3_J_2_15_P_2_30_kT_5.3_DGDK_0.05.mat';
% datafname = 'p2_linear/8x24_ky_0.3_J_2_15_P_2_20_kT_6.96_DGDK_0.02.mat';
% datafname = 'p2_linear/8x24_ky_0.3_J_2_15_P_2_30_kT_5.3_DGDK_0.025.mat';
%% new
% datafname = 'p2_linear_new/8x24_ky_0.3_P_4_J_2_kT_3.5_6.96_nu_0_0.5_IVDK.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_P_4_J_2_kT_3.5_6.96_nu_0_0.5_DGDK.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_P_4_J_2_kT_3.5_6.96_nu_0_0.5_DGGK.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_P_4_J_2_kT_3.5_6.96_nu_0_0.5_SGGK.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_P_4_J_2_kT_3.5_6.96_nu_0_0.5_LDGK.mat';

% datafname = 'p2_linear_new/8x24_ky_0.3_P_6_J_3_kT_3.5_6.96_nu_0_0.5_DGGK.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_P_6_J_3_kT_3.5_6.96_nu_0_0.5_SGGK.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_P_6_J_3_kT_3.5_6.96_nu_0_0.5_LDGK.mat';

% datafname = 'p2_linear_new/8x24_ky_0.3_P_8_J_4_kT_3.5_6.96_nu_0_0.5_DGGK.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_P_8_J_4_kT_3.5_6.96_nu_0_0.5_SGGK.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_P_8_J_4_kT_3.5_6.96_nu_0_0.5_LDGK.mat';

% datafname = 'p2_linear_new/8x24_ky_0.3_kT_3_6.96_P_2_10_DGDK_0.05.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_kT_3_6.96_P_2_30_DGGK_0.05.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_kT_3_6.96_P_2_20_SGGK_0.05.mat';
datafname = 'p2_linear_new/8x24_ky_0.3_kT_3_6.96_P_2_10_LDGK_0.05.mat';

% datafname = 'p2_linear_new/8x24_ky_0.3_kT_3_6.96_P_2_20_DGDK_0.1.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_kT_3_6.96_P_2_20_DGGK_0.1.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_kT_3_6.96_P_2_20_SGGK_0.1.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_kT_3_6.96_P_2_10_LDGK_0.1.mat';

% datafname = 'p2_linear_new/8x24_ky_0.3_P_2_J_1_kT_3_5_nu_0.001_0.1_DGGK.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_P_4_J_2_kT_2.5_4.5_nu_0.001_0.1_DGGK.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_P_6_J_3_kT_2.5_4.5_nu_0.001_0.1_DGGK.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_P_8_J_4_kT_2.5_4.5_nu_0.001_0.1_DGGK.mat';
% datafname = 'p2_linear_new/8x24_ky_0.3_P_16_J_8_kT_2.5_4.5_nu_0.001_0.1_DGGK.mat';

% datafname = 'p2_linear_new/8x24_ky_0.05_1_P_2_64_DGDK_0.001_kT_6.96.mat';
% datafname = 'p2_linear_new/8x24_ky_0.05_1_P_2_64_DGDK_0.001_kT_5.3.mat';
%% ky pj scan
% datafname = 'p2_linear_new/8x24_ky_0.05_1_P_2_64_DGDK_0.001_kT_6.96.mat';
% datafname = 'p2_linear_new/8x24_ky_0.05_1_P_2_48_DGDK_0.001_kT_6.96.mat';
% datafname = 'p2_linear_new/8x24_ky_0.05_1_P_2_16_DGGK_0.01_kT_5.3.mat';
% datafname = 'p2_linear_new/8x24_ky_0.05_1_P_2_16_SGGK_0.01_kT_5.3.mat';
% datafname = 'p2_linear_new/8x24_ky_0.05_1_P_2_10_LDGK_0.01_kT_5.3.mat';

% datafname = 'p2_linear_new/8x24_P_4_ky_0.05_1_DGGK_nu_0.005_0.05_kT_5.3.mat';
% datafname = 'p2_linear_new/8x24_P_8_ky_0.05_1_DGGK_nu_0.005_0.05_kT_5.3.mat';
% datafname = 'p2_linear_new/8x24_P_16_ky_0.05_1_DGGK_nu_0.005_0.05_kT_5.3.mat';

% datafname = 'p2_linear_new/8x24_P_4_ky_0.05_1_SGGK_nu_0.005_0.05_kT_5.3.mat';
% datafname = 'p2_linear_new/8x24_P_8_ky_0.05_1_SGGK_nu_0.005_0.05_kT_5.3.mat';
% datafname = 'p2_linear_new/8x24_P_16_ky_0.05_1_SGGK_nu_0.005_0.05_kT_5.3.mat';

% datafname = 'p2_linear_new/8x24_P_4_ky_0.05_1_LDGK_nu_0.005_0.05_kT_5.3.mat';
% datafname = 'p2_linear_new/8x24_P_8_ky_0.05_1_LDGK_nu_0.005_0.05_kT_5.3.mat';
% datafname = 'p2_linear_new/8x24_P_10_ky_0.05_1_LDGK_nu_0.005_0.05_kT_5.3.mat';
% datafname = 'p2_linear_new/8x24_P_16_ky_0.05_1_LDGK_nu_0.005_0.05_kT_5.3.mat';

%% Chose if we filter gamma>0.05
FILTERGAMMA = 1;

%% Load data
fname = ['../results/',datafname];
d = load(fname);
if FILTERGAMMA
    d.data = d.data.*(d.data>0.025);
    d.err  = d.err.*(d.data>0.025);
end
if 1
%% Pcolor of the peak
figure;
% [XX_,YY_] = meshgrid(d.s1,d.s2);
[XX_,YY_] = meshgrid(1:numel(d.s1),1:numel(d.s2));
% pclr=imagesc_custom(XX_,YY_,d.data'.*(d.data>0)');
pclr=contourf(1:numel(d.s1),1:numel(d.s2),d.data'.*(d.data>0)');
% pclr=surf(1:numel(d.s1),1:numel(d.s2),d.data'.*(d.data>0)');
title(d.title);
xlabel(d.s1name); ylabel(d.s2name);
set(gca,'XTick',1:numel(d.s1),'XTicklabel',d.s1)
set(gca,'YTick',1:numel(d.s2),'YTicklabel',d.s2)
colormap(jet)
colormap(bluewhitered)
clb=colorbar; 
clb.Label.String = '$\gamma c_s/R$';
clb.Label.Interpreter = 'latex';
clb.Label.FontSize= 18;
end
if 0
%% Scan along first dimension
figure
colors_ = jet(numel(d.s2));
for i = 1:numel(d.s2)
    % plot(d.s1,d.data(:,i),'s-',...
    plot(d.s1(d.data(:,i)>0),d.data((d.data(:,i)>0),i),'s-',...
        'LineWidth',2.0,...
        'DisplayName',[d.s2name,'=',num2str(d.s2(i))],...
        'color',colors_(i,:)); 
    % errorbar(d.s1,d.data(:,i),d.err(:,i),'s-',...
    % 'LineWidth',2.0,...
    % 'DisplayName',[d.s2name,'=',num2str(d.s2(i))],...
    % 'color',colors_(i,:)); 
    hold on;
end
xlabel(d.s1name); ylabel(d.dname);title(d.title);
xlim([d.s1(1) d.s1(end)]);
colormap(colors_);
clb = colorbar;
% caxis([d.s2(1)-0.5,d.s2(end)+0.5]);
clim([1 numel(d.s2)+1]);
clb.Ticks=linspace(d.s2(1),d.s2(end),numel(d.s2));
clb.Ticks    =1.5:numel(d.s2)+1.5;
clb.TickLabels=d.s2;
clb.Label.String = d.s2name;
clb.Label.Interpreter = 'latex';
clb.Label.FontSize= 18;
end
if 0
%% Scan along second dimension
figure
colors_ = jet(numel(d.s1));
for i = 1:numel(d.s1)
    plot(d.s2,d.data(i,:),'s-',...
        'LineWidth',2.0,...
        'DisplayName',[d.s1name,'=',num2str(d.s1(i))],...
        'color',colors_(i,:)); 
%     errorbar(d.s2,d.data(i,:),d.err(i,:),'s-',...
%         'LineWidth',2.0,...
%         'DisplayName',[d.s1name,'=',num2str(d.s1(i))],...
%         'color',colors_(i,:)); 
    hold on;
end
xlabel(d.s2name); ylabel(d.dname);title(d.title);
xlim([d.s2(1) d.s2(end)]);
colormap(jet(numel(d.s1)));
clb = colorbar;
caxis([d.s1(1)-0.5,d.s1(end)+0.5]);
clb.Ticks=linspace(d.s1(1),d.s1(end),numel(d.s1));
clb.YTick=d.s1;
clb.Label.String = d.s1name;
clb.TickLabelInterpreter = 'latex';
clb.Label.Interpreter = 'latex';
clb.Label.FontSize= 18;
end

if 0
%% Convergence analysis
figure
% target_ = 0.25*(d.data(end,end)+d.data(end-i_,end)+d.data(end,end-i_)+d.data(end-i_,end-i_));
colors_ = jet(numel(d.s1));
for i = 1:numel(d.s1)
%     target_ = d.data(i,end);
    target_ = 2.79666916212537142172e-01; % Value for nuDGDK = 0.05, kT=6.96, (40,20), Nkx=8
    if target_ > 0
    eps_    = abs(target_ - d.data(i,1:end-1))/abs(target_);
    semilogy(d.s2(1:end-1),eps_,'-s',...
        'LineWidth',2.0,...
        'DisplayName',[d.s1name,'=',num2str(d.s1(i))],...
        'color',colors_(i,:));
    hold on;
    end
end
xlabel(d.s2name); ylabel('$\epsilon_r$');title(d.title);
xlim([d.s2(1) d.s2(end)]);
colormap(colors_);
clb = colorbar;
caxis([d.s1(1)-0.5,d.s1(end)+0.5]);
clb.Ticks=linspace(d.s1(1),d.s1(end),numel(d.s1));
clb.YTick=d.s1;
clb.Label.String = d.s1name;
clb.TickLabelInterpreter = 'latex';
clb.Label.Interpreter = 'latex';
clb.Label.FontSize= 18;
grid on;
end
if 0
%% Pcolor of the error
figure;  i_ = 0;
% target_ = 2.72724991618068013377e-01; % Value for nuDGDK = 1.0, kT=6.96, (40,20), Nkx=8
% target_ = 2.79666916212537142172e-01; % Value for nuDGDK = 0.05, kT=6.96, (40,20), Nkx=8
% target_ = 2.71455e-01; % Value for nuDGDK = 0.001, kT=6.96, (40,20), Nkx=8
% target_ = 1.39427e-01; % Value for nuDGDK = 0.001, kT=5.3, (30,16), Nkx=8
target_ = 2.72510405826983714839e-01  % Value for nuDGDK = 0.001, kT=6.96, (50,25), Nkx=8
% target_ = 2.73048910051283844069e-01; % Value for nuDGDK = 0.0, kT=6.96, (40,20), Nkx=8
% target_ = 0.25*(d.data(end,end)+d.data(end-i_,end)+d.data(end,end-i_)+d.data(end-i_,end-i_));
% eps_    = log(abs(target_ - d.data)/abs(target_));
eps_    = max(-10,log(abs(target_ - d.data)/abs(target_)));
sign_   = 1;%sign(d.data - target_);
eps_ = d.data;
for i = 1:numel(d.s1)
    for j = 1:numel(d.s2)
        % target_ = d.data(i,end);
        % target_ = d.data(i,end);
        % target_ = d.data(end,j);
        eps_(i,j) = log(abs(target_ - d.data(i,j))/target);
        if target_ > 0
    %     eps_(i,:) = max(-12,log(abs(target_ - d.data(i,1:end))/abs(target_)));
    %     eps_(i,:) = log(abs(target_ - d.data(i,1:end))/abs(target_));
    %     eps_(i,:) = min(100,100*abs(target_ - d.data(i,1:end))/abs(target_));
        else
        end
    end
end
[XX_,YY_] = meshgrid(d.s1,d.s2);
[XX_,YY_] = meshgrid(1:numel(d.s1),1:numel(d.s2));
pclr=imagesc_custom(XX_,YY_,eps_'.*(d.data>0)'.*sign_');
% pclr=contourf(1:numel(d.s1),1:numel(d.s2),eps_'.*(d.data>0)'.*sign_',5);
% pclr=surf(1:numel(d.s1),1:numel(d.s2),eps_'.*(d.data>0)'.*sign_');
title(d.title);
xlabel(d.s1name); ylabel(d.s2name);
set(gca,'XTick',1:numel(d.s1),'XTicklabel',d.s1)
set(gca,'YTick',1:numel(d.s2),'YTicklabel',d.s2)
% colormap(jet)
colormap(bluewhitered)
% caxis([-10, 0]);
clb=colorbar; 
clb.Label.String = '$\log(\epsilon_r)$';
clb.Label.Interpreter = 'latex';
clb.Label.FontSize= 18;
end