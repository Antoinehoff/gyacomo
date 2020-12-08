default_plots_options

%% Strong scaling measurement

% Handwritten results for 512x256, P,J=2,1, Tmax = 5
Results_512_21.np    = [   1,    2,    4,    8,    12,   16,   20,   24];
Results_512_21.time  = [0000, 0000, 0000, 0000,  0000, 0000, 0000, 0000];

% Handwritten results for 512x256, P,J=3,2, Tmax = 2, mu=0, dt 1e-2?
Results_512_32.np    = [   1,    2,    4,    8,    12,   16,   20,   24];
Results_512_32.time  = [0000, 0000, 0000, 0000,  0000, 0000, 0000, 0000];

% Handwritten results for 1024x512, P,J=1,1, Tmax = 5 dt = 0.05, mu = 0
Results_1024_11.np    = [   1,    2,    4,    8,    12,   16,   20,   24];
Results_1024_11.time  = [0000, 0000, 0000, 0000,  0000, 0000, 0000, 0000];

% Handwritten results for 1024x512, P,J=2,2, Tmax = 2 dt = 0.05, mu = 0
Results_1024_22.np    = [   1,    2,    4,    8,    12,   16,   20,   24];
Results_1024_22.time  = [0000, 0000, 0000, 0000,  0000, 0000, 0000, 0000];

% Handwritten results for 1024x512, P,J=6,4, Tmax = 2 dt = 0.05, mu = 0
Results_1024_64.np    = [   1,    2,    4,    8,    12,   16,   20,   24];
Results_1024_64.time  = [0000, 0000, 0000, 0000,  0000, 0000, 0000, 0000];

%
fig = figure;

plot(1:24,1:24,'-k','DisplayName','Ideal')
hold on
% res = Results_256_21;
% plot(res.np,res.time(1)./(res.time),'o--','DisplayName','$256\times128$, $P,J=2,1$');
res = Results_512_21;
plot(res.np,res.time(1)./(res.time),'v-','DisplayName','$512\times256$, $P,J=2,1$');
res = Results_512_32;
plot(res.np,res.time(1)./(res.time),'>-','DisplayName','$512\times256$, $P,J=3,2$');
res = Results_1024_11;
plot(res.np,res.time(1)./(res.time),'o-','DisplayName','$1024\times512$, $P,J=1,1$');
res = Results_1024_22;
plot(res.np,res.time(1)./(res.time),'s-','DisplayName','$1024\times512$, $P,J=2,2$');xlim([1,max(res.np)]);
res = Results_1024_64;
plot(res.np,res.time(1)./(res.time),'d-','DisplayName','$1024\times512$, $P,J=6,4$');xlim([1,max(res.np)]);
xlabel('$N_p$'); ylabel('speedup')
xlim([1,24]); ylim([1,24])
legend('show')
title('Strong scaling')
grid on  
    

FIGNAME = '/home/ahoffman/HeLaZ/results/strong_scaling_new.pdf';
saveas(fig,FIGNAME);
disp(['Figure saved @ : ',FIGNAME])
