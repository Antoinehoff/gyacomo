default_plots_options

NPS = [01 02 04 08 12 16 20 24];
TIMES = NPS;
if 0
%% Load times
i_ = 1;
for np = NPS
    SIM_NAME = sprintf('Scaling_np%02d',np);
    hostfile = ['/marconi_scratch/userexternal/ahoffman/HeLaZ/results/',SIM_NAME,'/',BASIC.PARAMS];
    localfile= ['../results/',SIM_NAME,'/'];
    system(['scp -r ahoffman@login.marconi.cineca.it:',hostfile,' ',localfile]);
    filename = ['../results/',SIM_NAME,'/',BASIC.PARAMS,'/outputs_00.h5'];
    TIMES(i_)   = h5readatt(filename,'/data/input','cpu_time');
    i_ = i_ + 1;
end
TIMES
end
%% Strong scaling measurement (no diagnostics)

% Handwritten results for 512x256, P,J=2,1, Tmax = 10, mu=0, dt = 5e-2
Results_512_21.np    = NPS;
Results_512_21.time  = [799   422   206   106    81    72    66    82]; %tmax 10
% Results_512_21.time  = [0162, 0108, 0055, 0032,  0030, 0045, 0061, 0084]; %tmax 2

% Handwritten results for 512x256, P,J=3,2, Tmax = 10, mu=0, dt 5e-2
Results_512_32.np    = [   1,    2,    4,    8,    12,   16,   20,   24];
Results_512_32.time  = [2421  1241  0598  0309   0221  0188  0159  0148];

% Handwritten results for 1024x512, P,J=2,1, Tmax = 5 dt = 0.05, mu = 0
Results_1024_21.np    = [   1,    2,    4,    8,    12,   16,   20,   24];
% Results_1024_21.time  = [1920, 0000, 0563, 0306,  0247, 0240, 0000, 0000];
Results_1024_21.time  = [1920, 1260, 0556, 0289,  0219, 0189, 0172, 0167];

% Handwritten results for 1024x512, P,J=3,2, Tmax = 2 dt = 0.05, mu = 0
Results_1024_32.np    = [   1,    2,    4,    8,    12,   16,   20,   24];
Results_1024_32.time  = [2199, 1424, 0623, 0324,  0243, 0208, 0188, 0180];

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
res = Results_1024_21;
plot(res.np,res.time(1)./(res.time),'o-','DisplayName','$1024\times512$, $P,J=2,1$');
res = Results_1024_32;
plot(res.np,res.time(1)./(res.time),'s-','DisplayName','$1024\times512$, $P,J=3,2$');xlim([1,max(res.np)]);
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
