default_plots_options

% NPS = [01 02 04 08 12 16 20 24];
NPS = [02 04 10];
TIMES = NPS;
if 0
%% Load times
i_ = 1;
nn = 2;
for np = NPS
    SIM_NAME = sprintf('nn%02d_np%02d',nn,np)
    hostfile = ['/marconi_scratch/userexternal/ahoffman/HeLaZ/results/',SIM_NAME,'/',BASIC.PARAMS];
    localtarget= ['../results/',sprintf('Scaling/nn%02d_np%02d',nn,np),'/'];
    system(['scp -r ahoffman@login.marconi.cineca.it:',hostfile,' ',localtarget]);
    filename = ['../results/',sprintf('Scaling/nn%02d_np%02d',nn,np),'/',BASIC.PARAMS,'/outputs_00.h5'];
    TIMES(i_)   = h5readatt(filename,'/data/input','cpu_time');
    i_ = i_ + 1;
end
disp(PARAMS)
disp(TIMES')
end
%% Strong scaling measurement (no diagnostics)

% Handwritten results for 512x256, P,J=2,1, Tmax = 10, mu=0, dt = 5e-2
Results_512_21.np    = NPS;
Results_512_21.tproc = [857   447   217   107    74    59    46    41]; %Nthreads
Results_512_21.tnode = [854   442   199    96    75    64    57    55]; %Nnodes

% Handwritten results for 512x256, P,J=3,2, Tmax = 10, mu=0, dt 5e-2
Results_512_32.np    = [   1,    2,    4,    8,    12,   16,   20,   24];
Results_512_32.tproc = [2546  1290  0630  0316   0220  0177  0143  0125];

% Handwritten results for 512x256, P,J=6,4, Tmax = 5, mu=0, dt 5e-2
Results_512_64.np    = [   1,    2,    4,    8,    12,   16,   20,   24];
Results_512_64.tproc = [5801  2974  1467  0746   0507  0408  0320  0279];

% Handwritten results for 1024x512, P,J=2,1, Tmax = 5 dt = 0.05, mu = 0
Results_1024_21.np    = [   1,    2,    4,    8,    12,   16,   20,   24];
Results_1024_21.tproc = [2051, 1296, 0583, 0296,  0205, 0163, 0132, 0115];
Results_1024_21.tnode = [2052  1327  0511  0231   0157  0109  0099  0083]; %Nnodes

% Handwritten results for 1024x512, P,J=3,2, Tmax = 2 dt = 0.05, mu = 0
Results_1024_32.np    = [   1,    2,    4,    8,    12,   16,   20,   24];
Results_1024_32.tproc = [2248, 1526, 0644, 0326,  0233, 0182, 0148, 0130];

% Handwritten results for 1024x512, P,J=6,4, Tmax = 2 dt = 0.05, mu = 0
Results_1024_64.np    = [   1,    2,    4,    8,    12,   16,   20,   24];
Results_1024_64.tproc = [10330,6548, 2901, 1497,  1001, 0769, 0610, 0520];

%
fig = figure;

    subplot(121)
    plot(1:24,1:24,'-k','DisplayName','Ideal')
    hold on
    res = Results_512_21;
    plot(res.np,res.tproc(1)./(res.tproc),'v-','Color',line_colors(1,:),...
        'DisplayName','$512\times256$, $P,J=2,1$');

    res = Results_512_32;
    plot(res.np,res.tproc(1)./(res.tproc),'v-','Color',line_colors(2,:),...
        'DisplayName','$512\times256$, $P,J=3,2$');

    res = Results_512_64;
    plot(res.np,res.tproc(1)./(res.tproc),'v-','Color',line_colors(3,:),...
        'DisplayName','$512\times256$, $P,J=6,4$');

    res = Results_1024_21;
    plot(res.np,res.tproc(1)./(res.tproc),'^--','Color',line_colors(1,:),...
        'DisplayName','$1024\times512$, $P,J=2,1$');

    res = Results_1024_32;
    plot(res.np,res.tproc(1)./(res.tproc),'^--','Color',line_colors(2,:),...
        'DisplayName','$1024\times512$, $P,J=3,2$');

    res = Results_1024_64;
    plot(res.np,res.tproc(1)./(res.tproc),'^--','Color',line_colors(3,:),...
        'DisplayName','$1024\times512$, $P,J=6,4$');

    xlabel('$N_p$'); ylabel('speedup')
    xlim([1,24]);
    legend('show')
    title('$N_n=01$')
    
subplot(122)
    plot(1:24,1:24,'-k','DisplayName','Ideal')
    hold on
    res = Results_512_21;
    plot(res.np,res.tnode(1)./(res.tnode),'o-','Color',line_colors(1,:),...
        'DisplayName','$512\times256$, $P,J=2,1$');

    res = Results_1024_21;
    plot(res.np,res.tnode(1)./(res.tnode),'o--','Color',line_colors(1,:),...
        'DisplayName','$1024\times512$, $P,J=2,1$');

    xlabel('$N_n$'); ylabel('speedup')
    xlim([1,24]);
    legend('show')
    title('$N_p=01$')
        

FIGNAME = '/home/ahoffman/HeLaZ/results/strong_scaling_new.pdf';
saveas(fig,FIGNAME);
disp(['Figure saved @ : ',FIGNAME])
