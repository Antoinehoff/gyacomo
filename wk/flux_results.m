default_plots_options
if 0
%% Compute time average and std of the mean flow
t0 = 160; t1 = 200; [~,it0] = min(abs(t0-Ts2D)); [~,it1] = min(abs(t1-Ts2D)); 
range  = it0:it1;
avg    = mean(Flux_ri(range))
stdev  = std(Flux_ri(range))^(.5)
figure
hist(Flux_ri(range),20); xlabel('$\Gamma$')
end
%% Handwritten results
% High definition results (256x128)
Results_256x128.Gamma = [0.62, 0.60,  2.7, 0.08];
Results_256x128.N     = [ 256,  256,  256,  256];
Results_256x128.L     = [  50,   66,   66,   50];
Results_256x128.P     = [   2,    2,    2,    3];
Results_256x128.J     = [   1,    1,    1,    2];
Results_256x128.NU    = [0.01, 0.01, 0.01, 0.01];
Results_256x128.etaB  = [ 0.5,  0.5,  0.4,  0.5];
Results_256x128.mrkr  = [ 'v',  'v',  'v',  'v'];
Results_256x128.clr   = [ 'b',  'b',  'b',  'g'];
% Low definition results (128x64)
Results_128x64.Gamma = [0.29, 0.05, 7e-4, 0.31,  3.7,  2e-3];
Results_128x64.N     = [128,   128,  128,  128   128,  128];
Results_128x64.L     = [  25,   25,   25,   33    33,   33];
Results_128x64.P     = [   2,    2,    2,    2,    2,    2];
Results_128x64.J     = [   1,    1,    1,    1,    1,    1];
Results_128x64.NU    = [0.01,  0.1, 0.01, 0.01, 0.01, 0.01];
Results_128x64.etaB  = [ 0.5,  0.5, 0.67,  0.5   0.4,  0.6];
Results_128x64.mrkr  = [ 's',  's',  's',  's',  's',  's'];
Results_128x64.clr   = [ 'b',  'b',   'b',  'r',  'r',  'r'];

Ricci_Rogers.Gamma = [ 1 1e-2];
Ricci_Rogers.etaB  = [0.5 1.0];

if 1
% Fig 3 of Ricci Rogers 2006
fig = figure;

res = Results_256x128;
for i = 1:numel(res.Gamma)
    if res.NU(i) == 0.01
    semilogy(res.etaB(i),res.Gamma(i),...
        res.mrkr(i),'DisplayName','256x128', 'color', res.clr(i));
    end
    hold on;
end
res = Results_128x64;
for i = 1:numel(res.Gamma)
    if res.NU(i) == 0.01
    semilogy(res.etaB(i),res.Gamma(i),...
        res.mrkr(i),'DisplayName','128x64', 'color', res.clr(i)); 
    end
    hold on;
end
   xlabel('$\eta_B$'); ylabel('$\Gamma^\infty_{part}$') 
end
plot(Ricci_Rogers.etaB,Ricci_Rogers.Gamma,'--ok');
grid on; title('$\nu = 0.01$')

FIGNAME = [SIMDIR,'flux_study.png'];
saveas(fig,FIGNAME);
disp(['Figure saved @ : ',FIGNAME])