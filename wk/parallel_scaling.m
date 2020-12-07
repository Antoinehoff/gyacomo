default_plots_options

%% Strong scaling measurement

% Handwritten results for 256x128, P,J=2,1, Tmax = 20
Results_256_21.np    = [   1,    2,    4,    8,    10,    12];
Results_256_21.time  = [2450, 1346, 0680, 0389,  0323,  0307];

% Handwritten results for 512x256, P,J=2,1, Tmax = 5
Results_512_21.np    = [   1,    2,    4,    8,    16,   20,   24];
Results_512_21.time  = [3429, 1680, 0842, 0443,  0292, 0322, 0362];

% Handwritten results for 512x256, P,J=3,2, Tmax = 2
Results_512_32.np    = [   1,    2,    4,    8,    16,   20];
Results_512_32.time  = [4450, 2267, 1136, 0595,  0363, 0323];

% Handwritten results for 1024x512, P,J=1,1, Tmax = 5 dt = 0.05, mu = 0
Results_1024_11.np    = [   1,    2,    4,    6,    8,   12,   16];
Results_1024_11.time  = [1568, 1046, 0490, 0347, 0257, 0221,  219];

% Handwritten results for 1024x512, P,J=2,2, Tmax = 2 dt = 0.05, mu = 0
Results_1024_22.np    = [   1,    2,    4,    6,    8,   10,   12,   16,   20];
Results_1024_22.time  = [2391, 1373, 0654, 0457, 0343, 0297, 0274, 0219, 0206];

% Handwritten results for 1024x512, P,J=6,4, Tmax = 2 dt = 0.05, mu = 0
Results_1024_22.np    = [   1,    2,    4,    6,    8,   10,   12,   16,   20];
Results_1024_22.time  = [2391, 1373, 0654, 0457, 0343, 0297, 0274, 0219, 0206];

%
fig = figure;

plot(1:24,1:24,'--k','DisplayName','Ideal')
hold on
res = Results_256_21;
plot(res.np,res.time(1)./(res.time),'o-','DisplayName','$256\times128$, $P,J=2,1$');
res = Results_512_21;
plot(res.np,res.time(1)./(res.time),'o-','DisplayName','$512\times256$, $P,J=2,1$');
res = Results_512_32;
plot(res.np,res.time(1)./(res.time),'o-','DisplayName','$512\times256$, $P,J=3,2$');
res = Results_1024_11;
plot(res.np,res.time(1)./(res.time),'o-','DisplayName','$1024\times512$, $P,J=1,1$');
res = Results_1024_22;
plot(res.np,res.time(1)./(res.time),'o-','DisplayName','$1024\times512$, $P,J=2,2$');xlim([1,max(res.np)]);
xlabel('$N_p$'); ylabel('speedup') 
legend('show')
title('Strong scaling')
grid on  
    

FIGNAME = [SIMDIR,'strong_scaling.png'];
saveas(fig,FIGNAME);
disp(['Figure saved @ : ',FIGNAME])

%% Weak scaling
% Handwritten results for P,J=2,1, Tmax = 5, dt = 0.01, Nz = Nr
Results_1_64.np    = [   1,    2,    4,    8];
Results_1_64.Nr    = [  64,   90,  128,  180];
Results_1_64.time  = [0064, 0074, 0082, 0101];

% Handwritten results for P,J=2,1, Tmax = 5, dt = 0.01, Nz = 128
Results_1_128.np    = [   1,    2,    4,    8,   16];
Results_1_128.Nr    = [  32,   64,  128,  256,  512];
Results_1_128.time  = [0032, 0037, 0043, 0049, 0070];

% Handwritten results for Tmax = 5, dt = 0.05, mu = 0, etab =0, Pi=Ji=Pe=Je=1
Results_1_128.np    = [   1,    2,    4,    6,   16];
Results_1_128.N     = [ 256,  360,  512,  720, 1024];
Results_1_128.time  = [0059, 0072, 0000, 0153, 0070];

%
fig = figure;

plot(Results_1_64.np,Results_1_64.time,'ob','DisplayName','$256\times128$');
hold on
plot(Results_1_64.np,Results_1_64.time(1)*ones(numel(Results_1_64.np)),'--b','DisplayName','Ideal')

plot(Results_1_128.np,Results_1_128.time,'or','DisplayName','$256\times128$');
plot(Results_1_128.np,Results_1_128.time(1)*ones(numel(Results_1_128.np)),'--r','DisplayName','Ideal')

xlim([1,max(res.np)]);
xlabel('$N_p$'); ylabel('CPU time [s]') 
legend('show')
title('Weak scaling')
grid on  
    

FIGNAME = [SIMDIR,'weak_scaling.png'];
saveas(fig,FIGNAME);
disp(['Figure saved @ : ',FIGNAME])
