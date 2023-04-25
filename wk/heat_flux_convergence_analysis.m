if 0
    %% GYAC Local low res
    % resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/local_runs/nu=0.05/3x2x64x48x16/CBC/';
    % resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/local_runs/3x2x64x48x16/CBC/';
    % resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/local_runs/5x2x64x48x16/CBC/';
    resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/local_runs/7x3x64x48x16/CBC/';
    % resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/local_runs/9x2x64x48x16/CBC/';
    % resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/local_runs/11x2x64x48x16/CBC/';
    % resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/local_runs/17x2x64x48x16/CBC/';
    data = {};
    data    = compile_results_low_mem(data,resdir,00,10);
    % fast heat flux analysis
    T    = data.Ts0D;
    Qx   = data.HFLUX_X;
    [~,it0] = min(abs(0.25*T(end)-T));
    Qavg = mean(Qx(it0:end));
    Qstd = std(Qx(it0:end));
    figure
    plot(T,Qx); hold on
    plot([T(it0) T(end)],Qavg*[1 1],'--k','DisplayName',...
    ['$Q_{avg}=',sprintf('%2.2f',Qavg),'\pm',sprintf('%2.2f',Qstd),'$']);
    legend('show')
    disp(['$Q_{avg}=',sprintf('%2.2f',Qavg),'\pm',sprintf('%2.2f',Qstd),'$']);
end
if 1
%% Manually gathered data
QvsPJ = [...
   %vp mu  Qxav  Qxer   kT
    03 02 62.19 26.46 6.96;...
    05 02 51.49 08.38 6.96;...
    07 03 50.08 08.38 6.96;...
    09 02 42.49 08.52 6.96;...
    11 02 36.84 07.22 6.96;... % 192x96 instead of 128x64
    17 02 37.26 07.63 6.96;...
];
figure
errorbar(QvsPJ(:,1).*QvsPJ(:,2),QvsPJ(:,3),QvsPJ(:,4),'--s',...
    'DisplayName','GYAC 64x48x16 ($\nu_{DGDK}=0.001$)'); hold on
end    
if 0
%% Manually gathered data
QvsPJ = [... % Old results done with nu = 0.05
   %vp mu  Qxav  Qxer   kT
    03 02 97.38 17.42 6.96;...
    05 02 64.80 13.50 6.96;...
    09 02 54.19 10.56 6.96;...
    11 02 59.57 13.58 6.96;... % 192x96 instead of 128x64
];
errorbar(QvsPJ(:,1).*QvsPJ(:,2),QvsPJ(:,3),QvsPJ(:,4),'--s',...
    'DisplayName','GYAC 64x48x16 ($\nu_{DGDK}=0.05$)'); hold on
end    
%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    %% GYAC Marconi standard res   
% Marconi standard res    
    % resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/CBC/3x2x128x64x24/';
%     resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/CBC/3x2x128x64x24_nu_5e-2/';
%     resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/CBC/5x3x128x64x24/';
%     resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/CBC/7x4x128x64x24/';
%     resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/local_runs/7x4x128x64x24/CBC_L120/';
%     resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/local_runs/7x4x128x64x24/CBC_L180/';
%     resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/CBC/21x6x192x96x24/';
    data = {};
    data    = compile_results_low_mem(data,resdir,00,10);
    % fast heat flux analysis
    T    = data.Ts0D;
    Qx   = data.HFLUX_X;
    [~,it0] = min(abs(0.25*T(end)-T));
    Qavg = mean(Qx(it0:end));
    Qstd = std(Qx(it0:end));
    figure
    plot(T,Qx); hold on
    plot([T(it0) T(end)],Qavg*[1 1],'--k','DisplayName',...
    ['$Q_{avg}=',sprintf('%2.2f',Qavg),'\pm',sprintf('%2.2f',Qstd),'$']);
    legend('show')
    disp(['$Q_{avg}=',sprintf('%2.2f',Qavg),'\pm',sprintf('%2.2f',Qstd),'$']);
end
if 1
%% Manually gathered data
QvsPJ = [...
   %vp mu  Qxav  Qxer   kT
    03 02 43.60 11.47 6.96;...
    03 02 65.03 17.79 6.96;... % nuDGDK = 0.05 instead of 0.001
    05 03 44.08 06.51 6.96;...
    07 04 33.09 04.85 6.96;...
    07 04 35.31 04.15 6.96;... % other run
    07 04 37.60 04.65 6.96;... % L=180
    21 06 33.07 05.31 6.96;... % 192x96 instead of 128x64
    % 21 11 00.00 00.00 6.96;... % 192x96 instead of 128x64
];
%     figure
errorbar(QvsPJ(:,1).*QvsPJ(:,2),QvsPJ(:,3),QvsPJ(:,4),'--s',...
    'DisplayName','GYAC 128x64x24 ($\nu_{DGDK}=0.001$)')
end

%%%%%%%%%%%%%%%%%%%%%%
if 0
    %% GENE
    % folder = '/misc/gene_results/CBC/KT_6.96_64x32x32x24x12_Nexc_5/';
%     folder = '/misc/gene_results/CBC/KT_6.96_128x64x24x8x4_Nexc_5_00/';
%     folder = '/misc/gene_results/CBC/KT_6.96_128x64x24x16x8_Nexc_5_00/';
%     folder = '/misc/gene_results/CBC/KT_6.96_128x64x24x32x16_Nexc_5_00/';
%     folder = '/misc/gene_results/CBC/KT_6.96_128x64x24x32x16_Nexc_5_01/';

    % folder = '/misc/gene_results/CBC/KT_5.3_128x64x24x32x16_Nexc_5_00/';
    % folder = '/misc/gene_results/CBC/KT_5.3_128x64x24x32x16_Nexc_5_01/';
    % folder = '/misc/gene_results/CBC/new_sim/KT_5.3_128x64x24x16x8_Nexc_5/';
    % folder = '/misc/gene_results/CBC/new_sim/KT_5.3_128x64x24x8x4_Nexc_5/';
    folder = '/misc/gene_results/CBC/new_sim/KT_6.96_128x64x24x8x4_Nexc_5_smallvbox/';
    % folder = '/misc/gene_results/CBC/new_sim/KT_6.96_128x64x24x16x8_Nexc_5_largexbox/';
    % folder = '/misc/gene_results/CBC/KT_5.3_128x64x24x16x8_Muz_0.02/';


    % folder = '/misc/gene_results/CBC/128x64x16x6x4/';
    % fast heat flux analysis
    nrgfile           = 'nrg.dat.h5';
    % nrgfile           = 'nrg_1.h5';
    T    = h5read([folder,nrgfile],'/nrgions/time');
    Qx   = h5read([folder,nrgfile],'/nrgions/Q_es');
    [~,it0] = min(abs(0.25*T(end)-T));
    Qavg = mean(Qx(it0:end));
    Qstd = std(Qx(it0:end));
    figure
    plot(T,Qx); hold on
    plot([T(it0) T(end)],Qavg*[1 1],'--k','DisplayName',...
    ['$Q_{avg}=',sprintf('%2.2f',Qavg),'\pm',sprintf('%2.2f',Qstd),'$']);
    legend('show')
    disp(['$Q_{avg}=',sprintf('%2.2f',Qavg),'\pm',sprintf('%2.2f',Qstd),'$']);
end
if 1
%% Manually gathered data
    QvsNv = [...
       %vp mu  Qxav Qxer   kT
        06 04 17.90 6.85 6.96;...
        08 04 02.01 0.59 6.96;...
        08 04 13.12 3.20 6.96;...
        16 08 38.13 6.08 6.96;...
        16 08 35.74 4.85 6.96;...
        32 16 33.10 7.01 6.96;...
    ];
%         figure
    errorbar((QvsNv(:,1).*QvsNv(:,2)),QvsNv(:,3),QvsNv(:,4),'--sk',...
        'DisplayName','GENE 128x64x24')

%         QvsNv = [...
%            %vp mu  Qxav Qxer   kT
%             32 16 00.32 0.09 5.30;...
%             16 08 00.77 0.18 5.30;...
%         ];
%         figure
%         errorbar(QvsNv(:,1).*QvsNv(:,2)/2,QvsNv(:,3),QvsNv(:,4),'sk')
%         legend('GENE KT=5.3');
end
set(gca,'xscale','log');
xlabel('$N_{v\parallel}\times N_{\mu}$'); ylabel('$Q_x$');
top_title('collisionless CBC');