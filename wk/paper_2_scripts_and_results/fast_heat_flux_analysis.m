%% Analysis of sequential kT scans
% resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/collisionless/kT_scan_nu_1e-3/5x3x128x64x24_dp/';
% resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/collisionless/kT_scan_nu_1e-3/7x4x128x64x24_dp/';
resdir = '/misc/gyacomo23_outputs/paper_2_GYAC23/collisionless/kT_scan_nu_1e-3/9x5x128x64x24_dp/';
Jobs    = [0 2 3 4 5];
kTA     = 0*Jobs;
kNA     = 0*Jobs;
QxstdA  = kTA;
QxavgA  = kTA;
    
figure

for ij = 1:numel(Jobs)
    data = {};
    data    = compile_results_low_mem(data,resdir,Jobs(ij),Jobs(ij));
    % fast heat flux analysis
    T    = data.Ts0D-data.Ts0D(1);
    Qx   = data.HFLUX_X;
    [~,it0] = min(abs(0.25*T(end)-T));
    Qavg = mean(Qx(it0:end));
    Qstd = std(Qx(it0:end));
    kTA(ij)    = data.inputs.K_T;
    kNA(ij)    = data.inputs.K_N;
    QxavgA(ij) = Qavg;
    QxstdA(ij) = Qstd;
    subplot(121)
    plot(T,Qx,'DisplayName',['$\kappa_T=',num2str(data.inputs.K_T),'$']); hold on
    plot([T(it0) T(end)],Qavg*[1 1],'--k','DisplayName',...
    ['$Q_{avg}=',sprintf('%2.2f',Qavg),'\pm',sprintf('%2.2f',Qstd),'$']);
end
[~,indices] = sort(kTA);
kTA = kTA(indices);
kNA = kNA(indices);
QxavgA = QxavgA(indices);
QxstdA = QxstdA(indices);

% Plots

legend('show')
subplot(122)
errorbar(kTA,QxavgA./kTA./kNA,QxstdA./kTA./kNA,'--sr',...
    'DisplayName','GYAC LD CBC ($\nu_{DGDK}=0.05$)'); hold on
xlabel('$\kappa_T$'); ylabel('$\chi$')
