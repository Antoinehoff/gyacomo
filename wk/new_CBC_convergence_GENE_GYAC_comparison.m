figure 
hold on  
TRTW = [0.5 1];
ERRBAR = 1;
%%%%%%%%%%%%%%%%%%%%%%%%
%% GYAC Marconi standard res   
resdirs = {...
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_1e-3/3x2x128x64x24/kT_7.0/';
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_1e-3/5x3x128x64x24/kT_7.0/';
     '/misc/gyacomo23_outputs/paper_2_GYAC23/save/kT_scan_1e-3/7x4x128x64x24/kT_7.0/';
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_1e-3/9x5x128x64x24/kT_7.0/';
     '/misc/gyacomo23_outputs/paper_2_GYAC23/save/kT_scan_1e-3/11x6x128x64x24/kT_7.0/';
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_1e-3/17x9x128x64x24/kT_7.0/';
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_1e-3/31x16x128x64x24/kT_7.0/'
    };
QvsPJ = zeros(numel(resdirs),4);
N = numel(resdirs);
for i = 1:N
    data = {};
    data    = compile_results_low_mem(data,resdirs{i},00,10);
    % fast heat flux analysis
    T    = data.Ts0D;
    Qx   = data.HFLUX_X;
    Trange = T(end)*TRTW;
    [~,it0] = min(abs(Trange(1)-T));
    [~,it1] = min(abs(Trange(2)-T));
    Qavg = mean(Qx(it0:it1));
    Qstd = std(Qx(it0:it1));
    disp(['$Q_{avg}=',sprintf('%2.2f',Qavg),'\pm',sprintf('%2.2f',Qstd),'$']);
    QvsPJ(i,1) = data.grids.Np;
    QvsPJ(i,2) = data.grids.Nj;
    QvsPJ(i,3) = Qavg;
    QvsPJ(i,4) = Qstd;
end
if ERRBAR
    errorbar(QvsPJ(:,1).*QvsPJ(:,2),QvsPJ(:,3),QvsPJ(:,4),'LineStyle','-',...
    'DisplayName','GM, $\nu_{DGDK}=0.001$','Marker','o',...
            'MarkerSize',7,'LineWidth',2);
else
    plot((QvsPJ(:,1).*QvsPJ(:,2)),QvsPJ(:,3),'-o','DisplayName','GM, $\nu_{DGDK}=0.001$')
end

%% GYAC Marconi nu = 0.01
resdirs = {...
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_1e-2/5x3x128x64x24/kT_7.0/';
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_1e-2/7x4x128x64x24/kT_7.0/';
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_1e-2/9x5x128x64x24/kT_7.0/';
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_1e-2/11x6x128x64x24/kT_7.0/';
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_1e-2/13x7x128x64x24/kT_7.0/'
    };
QvsPJ = zeros(numel(resdirs),4);
N = numel(resdirs);
for i = 1:N
    data = {};
    data    = compile_results_low_mem(data,resdirs{i},00,10);
    % fast heat flux analysis
    T    = data.Ts0D;
    Qx   = data.HFLUX_X;
    Trange = T(end)*TRTW;
    [~,it0] = min(abs(Trange(1)-T));
    [~,it1] = min(abs(Trange(2)-T));
    Qavg = mean(Qx(it0:it1));
    Qstd = std(Qx(it0:it1));
    disp(['$Q_{avg}=',sprintf('%2.2f',Qavg),'\pm',sprintf('%2.2f',Qstd),'$']);
    QvsPJ(i,1) = data.grids.Np;
    QvsPJ(i,2) = data.grids.Nj;
    QvsPJ(i,3) = Qavg;
    QvsPJ(i,4) = Qstd;
end
if ERRBAR
    errorbar(QvsPJ(:,1).*QvsPJ(:,2),QvsPJ(:,3),QvsPJ(:,4),'LineStyle','none',...
    'DisplayName','GM, $\nu_{DGDK}=0.01$','Marker','o',...
            'MarkerSize',7,'LineWidth',2);
else
    plot((QvsPJ(:,1).*QvsPJ(:,2)),QvsPJ(:,3),'-o','DisplayName','GM, $\nu_{DGDK}=0.01$')
end

%% GYAC Marconi nu = 0.05
resdirs = {...
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_5e-2/3x2x128x64x24/kT_7.0/';
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_5e-2/5x3x128x64x24/kT_7.0/';
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_5e-2/7x4x128x64x24/kT_7.0/';
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_5e-2/9x5x128x64x24/kT_7.0/';
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_5e-2/11x6x128x64x24/kT_7.0/';
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_5e-2/13x7x128x64x24/kT_7.0/'
     '/misc/gyacomo23_outputs/paper_2_GYAC23/kT_scan_nu_5e-2/15x8x128x64x24/kT_7.0/'
    };
QvsPJ = zeros(numel(resdirs),4);
N = numel(resdirs);
for i = 1:N
    data = {};
    data    = compile_results_low_mem(data,resdirs{i},00,10);
    % fast heat flux analysis
    T    = data.Ts0D;
    Qx   = data.HFLUX_X;
    Trange = T(end)*TRTW;
    [~,it0] = min(abs(Trange(1)-T));
    [~,it1] = min(abs(Trange(2)-T));
    Qavg = mean(Qx(it0:it1));
    Qstd = std(Qx(it0:it1));
    disp(['$Q_{avg}=',sprintf('%2.2f',Qavg),'\pm',sprintf('%2.2f',Qstd),'$']);
    QvsPJ(i,1) = data.grids.Np;
    QvsPJ(i,2) = data.grids.Nj;
    QvsPJ(i,3) = Qavg;
    QvsPJ(i,4) = Qstd;
end
if ERRBAR
    errorbar(QvsPJ(:,1).*QvsPJ(:,2),QvsPJ(:,3),QvsPJ(:,4),'LineStyle','none',...
    'DisplayName','GM, $\nu_{DGDK}=0.05$','Marker','o',...
            'MarkerSize',7,'LineWidth',2);
else
    plot((QvsPJ(:,1).*QvsPJ(:,2)),QvsPJ(:,3),'-o','DisplayName','GM, $\nu_{DGDK}=0.05$')
end
%% GENE normal box
resdirs = {...
     '/misc/gene_results/kT_scan_nu0/8x4x128x64x24/kT_7.0/';
     '/misc/gene_results/kT_scan_nu0/16x8x128x64x24/kT_7.0/';
     '/misc/gene_results/kT_scan_nu0/30x16x128x64x24/kT_7.0/';
     '/misc/gene_results/kT_scan_nu0/42x24x128x64x24/kT_7.0/'
    };
QvsNv = zeros(numel(resdirs),4);
N = numel(resdirs);
for i = 1:N
    % fast heat flux analysis
    nrgfile           = 'nrg.dat.h5';
    % nrgfile           = 'nrg_1.h5';
    T    = h5read([resdirs{i},nrgfile],'/nrgions/time');
    Qx   = h5read([resdirs{i},nrgfile],'/nrgions/Q_es');
    Trange = T(end)*TRTW;
    [~,it0] = min(abs(Trange(1)-T));
    [~,it1] = min(abs(Trange(2)-T));
    Qavg = mean(Qx(it0:end));
    Qstd = std(Qx(it0:end));
    disp(['$Q_{avg}=',sprintf('%2.2f',Qavg),'\pm',sprintf('%2.2f',Qstd),'$']);
    namelist      = read_namelist([resdirs{i},'parameters']);
    QvsNv(i,1) = namelist.box.nv0;
    QvsNv(i,2) = namelist.box.nw0;
    QvsNv(i,3) = Qavg;
    QvsNv(i,4) = Qstd;
end
if ERRBAR
    errorbar((QvsNv(:,1).*QvsNv(:,2)),QvsNv(:,3),QvsNv(:,4),'k','LineStyle','-',...
        'DisplayName','GENE','Marker','s',...
            'MarkerSize',10,'LineWidth',2);
else
    plot((QvsNv(:,1).*QvsNv(:,2)),QvsNv(:,3),'sk','DisplayName','GENE')
end
%% GENE half vbox
resdirs = {...
     '/misc/gene_results/kT_scan_nu0/8x4x128x64x24_half_vbox/kT_7.0/';
     '/misc/gene_results/kT_scan_nu0/16x8x128x64x24_half_vbox/kT_7.0/';
    };
QvsNv = zeros(numel(resdirs),4);
N = numel(resdirs);
for i = 1:N
    % fast heat flux analysis
    nrgfile           = 'nrg.dat.h5';
    % nrgfile           = 'nrg_1.h5';
    T    = h5read([resdirs{i},nrgfile],'/nrgions/time');
    Qx   = h5read([resdirs{i},nrgfile],'/nrgions/Q_es');
    [~,it0] = min(abs(0.25*T(end)-T));
    Qavg = mean(Qx(it0:end));
    Qstd = std(Qx(it0:end));
    disp(['$Q_{avg}=',sprintf('%2.2f',Qavg),'\pm',sprintf('%2.2f',Qstd),'$']);
    namelist      = read_namelist([resdirs{i},'parameters']);
    QvsNv(i,1) = namelist.box.nv0;
    QvsNv(i,2) = namelist.box.nw0;
    QvsNv(i,3) = Qavg;
    QvsNv(i,4) = Qstd;
end
if ERRBAR
    errorbar((QvsNv(:,1).*QvsNv(:,2)),QvsNv(:,3),QvsNv(:,4),'k','LineStyle','none',...
        'DisplayName','GENE, $L_{v_\parallel}/2,L_\mu/2$','Marker','x',...
            'MarkerSize',7,'LineWidth',2);
else
    plot((QvsNv(:,1).*QvsNv(:,2)),QvsNv(:,3),'xk','DisplayName','GENE, $L_{v_\parallel}/2,L_\mu/2$')
end
%% plot formatting
set(gca,'xscale','log');
xlabel('$N_{v\parallel}\times N_{\mu}$'); ylabel('$Q_x$');
grid on