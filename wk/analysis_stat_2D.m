if 0
%% ni ne phase space for different time windows
fig = figure;
    t0 = 3500; t1 = 4500; [~,it0] = min(abs(t0-Ts2D)); [~,it1] = min(abs(t1-Ts2D)); 
    plt = @(x) abs(squeeze(reshape(x(2,:,it0:it1),[],1)));
    plot(plt(Ni00),plt(Ne00),'.');
    title(['$',num2str(Ts2D(it0)),'\leq t \leq',num2str(Ts2D(it1)),'$']);
    xlabel('$n_i$'); ylabel('$n_e$');
end

if 0
%% histograms
t0  = 2700; t1 = 2800; [~,it0] = min(abs(t0-Ts2D)); [~,it1] = min(abs(t1-Ts2D)); 
plt = @(x) squeeze(reshape(x(:,:,it0:it1),[],1));
fig = figure('Position',  [100, 100, 800, 500]);
subplot(221)
    hist(plt(ne00),50); hold on
    title(['$',num2str(Ts2D(it0)),'\leq t \leq',num2str(Ts2D(it1)),'$']);
    xlabel('$n_e$');
subplot(222)
    hist(plt(ni00),50); hold on
    title(['$',num2str(Ts2D(it0)),'\leq t \leq',num2str(Ts2D(it1)),'$']);
    xlabel('$n_i$');
subplot(223)
    hist(plt(phi),50); hold on
    title(['$',num2str(Ts2D(it0)),'\leq t \leq',num2str(Ts2D(it1)),'$']);
    xlabel('$\phi$');
FIGNAME = ['histograms_ne_ni_phi'];
save_figure
end

if 0
%% ni ne phi 3D phase space
fig = figure;
    t0 = 2700; [~,it0] = min(abs(t0-Ts2D));
    plt = @(x) squeeze(reshape(x(:,:,it0),[],1));
    plot3(plt(ni00),plt(ne00),plt(phi),'.','MarkerSize',1.9);
    title(['$',num2str(Ts2D(it0)),'\leq t \leq',num2str(Ts2D(it1)),'$']);
    xlabel('$n_i$'); ylabel('$n_e$'); zlabel('$\phi$'); grid on;
end

t0    = 2400;
skip_ = 2; 
DELAY = 0.01*skip_;
FRAMES = floor(t0/(Ts2D(2)-Ts2D(1)))+1:skip_:numel(Ts2D);
if 0
%% ni ne phi 3D phase space GIF
GIFNAME = ['ni_ne_phi_phase_space',sprintf('_%.2d',JOBNUM)];
X = ne00; Y = ni00; Z = phi; T = Ts2D;
MARKERSIZE = 0.01;
XNAME = '$n_i$'; YNAME = '$n_e$'; ZNAME = '$\phi$'; VIEW = [1,-1,1];
create_gif_3D
end