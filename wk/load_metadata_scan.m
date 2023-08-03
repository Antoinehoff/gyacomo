% Metadata path
gyacomodir  = pwd;
gyacomodir = gyacomodir(1:end-2);
addpath(genpath([gyacomodir,'matlab'])) % ... add
addpath(genpath([gyacomodir,'matlab/plot'])) % ... add
addpath(genpath([gyacomodir,'matlab/compute'])) % ... add
addpath(genpath([gyacomodir,'matlab/load'])) % ... add% EXECNAME = 'gyacomo_1.0';

% datafname = 'lin_KBM/12x24_ky_0.05_0.75_P_2_16_DGGK_0.01_be_0.03.mat';
% datafname = 'lin_KBM/12x24_ky_0.05_0.75_P_2_16_DGGK_0.01_be_0.03.mat';
% datafname = 'lin_AUG_scan/12x24_ky_0.05_0.75_P_2_16_DGGK_0.01_be_0.000152.mat';
% datafname = 'lin_Entropy_scan/2x1_ky_0.05_0.75_P_2_8_DGDK_0_be_0.mat';
% datafname = 'lin_DTT_AB_rho85_PT_scan/16x24_ky_0.1_0.5_P_2_8_DGGK_0.05_be_0.0034.mat';
% datafname = 'lin_DTT_AB_rho85_PT_scan/16x24_ky_0.1_0.5_P_2_8_DGGK_0.05_be_0.0039.mat';
% datafname = 'lin_DTT_AB_rho85_PT_scan/16x24_ky_0.1_1_P_2_8_DGGK_0.05_be_0.0039.mat';
% datafname = 'lin_DTT_AB_rho85_PT_scan/16x24_ky_0.1_0.9_P_2_8_kN_1.33_DGGK_0.56901_be_0.0039.mat';
% datafname = 'lin_DTT_AB_rho85_PT_scan/8x24_ky_0.05_1.5_P_2_8_kN_1.33_DGGK_0.1_be_0.0039.mat';
% datafname = 'lin_JET_rho97_scan/4x32_ky_0.01_10_P_2_8_kN_10_DGGK_0.1_be_0.0031.mat';
% datafname = 'lin_JET_rho97_scan/4x32_ky_0.01_10_P_2_8_kN_10_DGGK_0.05_be_0.0031.mat';
% datafname = 'lin_JET_rho97_scan/8x32_ky_0.01_10_P_2_8_kN_10_DGGK_0.05_be_0.0031.mat';
% datafname = 'lin_JET_rho97_scan/8x32_ky_0.01_10_P_2_2_kN_10_DGGK_0.05_be_0.0031.mat';
datafname = 'lin_JET_rho97_scan/8x32_ky_0.031623_31.6228_P_2_2_kN_10_DGGK_0.1_be_0.0031.mat';
%% Chose if we filter gamma>0.05
FILTERGAMMA = 0;

%% Load data
fname = ['../results/',datafname];
d = load(fname);
gamma = real(d.data); g_err = real(d.err);
omega = imag(d.data); w_err = imag(d.err);
if FILTERGAMMA
    d.data = real(gamma).*(real(gamma)>0.025) + imag(gamma);
end
if 0
%% Pcolor of the peak
figure;
% [XX_,YY_] = meshgrid(d.s1,d.s2);
[XX_,YY_] = meshgrid(1:numel(d.s1),1:numel(d.s2));
pclr=imagesc_custom(XX_,YY_,gamma'.*(gamma>0)');
% pclr=contourf(1:numel(d.s1),1:numel(d.s2),gamma'.*(gamma>0)');
% pclr=surf(1:numel(d.s1),1:numel(d.s2),gamma'.*(gamma>0)');
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
if 1
%% Scan along first dimension
figure
colors_ = jet(numel(d.s2));
subplot(121)
for i = 1:numel(d.s2)
    % plot(d.s1,gamma(:,i),'s-',...
    % plot(d.s1(gamma(:,i)>0),gamma((gamma(:,i)>0),i),'s-',...
    %     'LineWidth',2.0,...
    %     'DisplayName',[d.s2name,'=',num2str(d.s2(i))],...
    %     'color',colors_(i,:)); 
    errorbar(d.s1,gamma(:,i),g_err(:,i),'s-',...
    'LineWidth',2.0,...
    'DisplayName',[d.s2name,'=',num2str(d.s2(i))],...
    'color',colors_(i,:)); 
    hold on;
end
xlabel(d.s1name); ylabel(d.dname);title(d.title);
xlim([d.s1(1) d.s1(end)]);

subplot(122)
for i = 1:numel(d.s2)
    % plot(d.s1,gamma(:,i),'s-',...
    % plot(d.s1(gamma(:,i)>0),gamma((gamma(:,i)>0),i),'s-',...
    %     'LineWidth',2.0,...
    %     'DisplayName',[d.s2name,'=',num2str(d.s2(i))],...
    %     'color',colors_(i,:)); 
    errorbar(d.s1,omega(:,i),w_err(:,i),'s-',...
    'LineWidth',2.0,...
    'DisplayName',[d.s2name,'=',num2str(d.s2(i))],...
    'color',colors_(i,:)); 
    hold on;
end
xlabel(d.s1name); ylabel('$\omega c_s/R$');title(d.title);
xlim([d.s1(1) d.s1(end)]);

colormap(colors_);
clb = colorbar;
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
    plot(d.s2,gamma(i,:),'s-',...
        'LineWidth',2.0,...
        'DisplayName',[d.s1name,'=',num2str(d.s1(i))],...
        'color',colors_(i,:)); 
%     errorbar(d.s2,gamma(i,:),g_err(i,:),'s-',...
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
% target_ = 0.25*(gamma(end,end)+gamma(end-i_,end)+gamma(end,end-i_)+gamma(end-i_,end-i_));
colors_ = jet(numel(d.s1));
for i = 1:numel(d.s1)
%     target_ = gamma(i,end);
    target_ = 2.79666916212537142172e-01; % Value for nuDGDK = 0.05, kT=6.96, (40,20), Nkx=8
    if target_ > 0
    eps_    = abs(target_ - gamma(i,1:end-1))/abs(target_);
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
% target_ = 0.25*(gamma(end,end)+gamma(end-i_,end)+gamma(end,end-i_)+gamma(end-i_,end-i_));
% eps_    = log(abs(target_ - gamma)/abs(target_));
eps_    = max(-10,log(abs(target_ - gamma)/abs(target_)));
sign_   = 1;%sign(gamma - target_);
eps_ = gamma;
for i = 1:numel(d.s1)
    for j = 1:numel(d.s2)
        % target_ = gamma(i,end);
        % target_ = gamma(i,end);
        % target_ = gamma(end,j);
        eps_(i,j) = log(abs(target_ - gamma(i,j))/target);
        if target_ > 0
    %     eps_(i,:) = max(-12,log(abs(target_ - gamma(i,1:end))/abs(target_)));
    %     eps_(i,:) = log(abs(target_ - gamma(i,1:end))/abs(target_));
    %     eps_(i,:) = min(100,100*abs(target_ - gamma(i,1:end))/abs(target_));
        else
        end
    end
end
[XX_,YY_] = meshgrid(d.s1,d.s2);
[XX_,YY_] = meshgrid(1:numel(d.s1),1:numel(d.s2));
pclr=imagesc_custom(XX_,YY_,eps_'.*(gamma>0)'.*sign_');
% pclr=contourf(1:numel(d.s1),1:numel(d.s2),eps_'.*(gamma>0)'.*sign_',5);
% pclr=surf(1:numel(d.s1),1:numel(d.s2),eps_'.*(gamma>0)'.*sign_');
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