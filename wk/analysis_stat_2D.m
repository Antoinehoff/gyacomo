if 0
%% ni ne phase space for different time windows
fig = figure;
    t0 = 20; t1 = 20; [~,it0] = min(abs(t0-Ts2D)); [~,it1] = min(abs(t1-Ts2D)); 
    plt = @(x) squeeze(reshape(x(:,:,it0:it1),[],1));
    plot(plt(ni00),plt(ne00),'.');
    title(['$',num2str(Ts2D(it0)),'\leq t \leq',num2str(Ts2D(it1)),'$']);
    xlabel('$n_i$'); ylabel('$n_e$');
end

if 0
%% histograms
fig = figure;
    t0 = 100; t1 = 100; [~,it0] = min(abs(t0-Ts2D)); [~,it1] = min(abs(t1-Ts2D)); 
    plt = @(x) squeeze(reshape(x(:,:,it0:it1),[],1));
    hist(plt(ni00),100); hold on
    title(['$',num2str(Ts2D(it0)),'\leq t \leq',num2str(Ts2D(it1)),'$']);
    xlabel('$n_i$');
end

if 0
%% ni ne phi 3D phase space
fig = figure;
    t0 = 20; t1 = 20; [~,it0] = min(abs(t0-Ts2D)); [~,it1] = min(abs(t1-Ts2D)); 
    plt = @(x) squeeze(reshape(x(:,:,it0:it1),[],1));
    plot3(plt(ni00),plt(ne00),plt(phi),'.','MarkerSize',0.5);
    title(['$',num2str(Ts2D(it0)),'\leq t \leq',num2str(Ts2D(it1)),'$']);
    xlabel('$n_i$'); ylabel('$n_e$'); zlabel('$\phi$'); grid on;
end

t0    = 0;
skip_ = 1; 
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