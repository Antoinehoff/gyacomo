default_plots_options
if 1
%% Compute time average and std of the mean flow
t0 = 180; t1 = 400; [~,it0] = min(abs(t0-Ts2D)); [~,it1] = min(abs(t1-Ts2D)); 
range  = it0:it1;
avg    = mean(Flux_ri(range))
stdev  = std(Flux_ri(range))^(.5)
figure
hist(Flux_ri(range),20); xlabel('$\Gamma$')
end
if 0
%% Handwritten results for nu = 0.01
% High definition results (256x128)
Results_256x128.Gamma = [0.02, 0.03, 0.20, 0.037,  2.7, 2.25,    4, 5e-4, 2e-3, 0.03];
Results_256x128.L     = [  66,   66,   66,    50,   66,   66,   66,   66,   66,   66];
Results_256x128.P     = [   2,    3,    4,     5,    2,    3,    4,    2,    3,    4];
Results_256x128.J     = [   1,    2,    2,     3,    1,    2,    2,    1,    2,    2];
Results_256x128.etaB  = [ 0.5,  0.5,  0.5,   0.5,  0.4,  0.4,  0.4,  0.6,  0.6,  0.6];
Results_256x128.mrkr  = [ 'v',  '>',  '^',   'o',  'v',  '>',  '^',  'v',  '>',  '^'];
Results_256x128.clr   = [ 'k',  'k',  'k',   'r',  'r',  'r',  'r',  'k',  'k',  'k'];
% Low definition results (128x64)
% Results_128x64.Gamma = [0.29, 0.05, 7e-4, 0.31,  3.7,  2e-3];
% Results_128x64.L     = [  25,   25,   25,   33    33,   33];
% Results_128x64.P     = [   2,    2,    2,    2,    2,    2];
% Results_128x64.J     = [   1,    1,    1,    1,    1,    1];
% Results_128x64.NU    = [0.01,  0.1, 0.01, 0.01, 0.01, 0.01];
% Results_128x64.etaB  = [ 0.5,  0.5, 0.67,  0.5   0.4,  0.6];
% Results_128x64.mrkr  = [ 's',  's',  's',  's',  's',  's'];
% Results_128x64.clr   = [ 'b',  'b',   'b',  'r',  'r',  'r'];

% Ricci_Rogers.Gamma = [2.5 1 1e-2];
% Ricci_Rogers.etaB  = [0.4 0.5 1.0];
Ricci_Rogers.Gamma = [10  1e-1];
Ricci_Rogers.etaB  = [0.5  1.0];
if 1
% Fig 3 of Ricci Rogers 2006
fig = figure;
semilogy(Ricci_Rogers.etaB,Ricci_Rogers.Gamma,'--','color',[0,0,0]+0.6);
hold on;
res = Results_256x128;
for i = 1:numel(res.Gamma)
    semilogy(res.etaB(i),res.Gamma(i),...
        res.mrkr(i),'DisplayName','256x128', 'color', res.clr(i));
    hold on;
end
% res = Results_128x64;
% for i = 1:numel(res.Gamma)
%     if res.NU(i) == 0.01
%     semilogy(res.etaB(i),res.Gamma(i),...
%         res.mrkr(i),'DisplayName','128x64', 'color', res.clr(i)); 
%     end
%     hold on;
% end
   xlabel('$\eta_B$'); ylabel('$\Gamma^\infty_{part}$') 
end
grid on; title('$\nu = 0.01$')
legend('Mix. Length, Ricci 2006','$P=2$, $J=1$','$P=3$, $J=2$','$P=4$, $J=2$','$P=5$, $J=3$')
FIGNAME = [SIMDIR,'flux_study_nu_1e-2.png']; ylim([1e-4, 10])
saveas(fig,FIGNAME);
disp(['Figure saved @ : ',FIGNAME])

%% Handwritten results for nu = 0.1
Results_256x128.Gamma = [0.026,0.026, 1e-2,    1,    1,    1, 2e-2,    1, 0.15,    3e-3];
Results_256x128.P     = [   2,     3,    4,    2,    3,    4,    2,    3,    4,       4];
Results_256x128.J     = [   1,     2,    2,    1,    2,    2,    1,    2,    2,       2];
Results_256x128.etaB  = [ 0.5,   0.5,  0.5,  0.4,  0.4,  0.4,  0.6,  0.6,  0.6,     0.7];
Results_256x128.mrkr  = [ 'v',   '>',  '^',  'v',  '>',  '^',  'v',  '>',  '^',     '^'];
Results_256x128.clr   = [ 'k',   'k',  'k',  'b',  'b',  'b',  'r',  'b',  'r',     'r'];

% Ricci_Rogers.Gamma = [2 1e-1];
% Ricci_Rogers.etaB  = [0.5 1.0];
Ricci_Rogers.Gamma = [10  1e-1];
Ricci_Rogers.etaB  = [0.5  1.0];

if 1
% Fig 3 of Ricci Rogers 2006
fig = figure;
semilogy(Ricci_Rogers.etaB,Ricci_Rogers.Gamma,'--','color',[0,0,0]+0.6);
hold on;
res = Results_256x128;
for i = 1:numel(res.Gamma)
    semilogy(res.etaB(i),res.Gamma(i),...
        res.mrkr(i),'DisplayName','256x128', 'color', res.clr(i));
    hold on;
end

   xlabel('$\eta_B$'); ylabel('$\Gamma^\infty_{part}$') 
end
grid on; title('$\nu = 0.1$')
legend('Mix. Length, Ricci 2006','$P=2$, $J=1$','$P=3$, $J=2$','$P=4$, $J=2$')
FIGNAME = [SIMDIR,'flux_study_nu_1e-2.png'];
saveas(fig,FIGNAME);
disp(['Figure saved @ : ',FIGNAME])
end