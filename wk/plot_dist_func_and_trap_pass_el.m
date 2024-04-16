[data.Napjz, data.Ts3D] = compile_results_3Da(DATADIR,J0,J1,'Napjz');
options.SHOW_FLUXSURF = 0;
options.SHOW_METRICS  = 0;
options.SHOW_CURVOP   = 0;
[fig, geo_arrays] = plot_metric(data,options);

eps_eff = 1/max((geo_arrays.hatB));
T = [40];
s = linspace(-4,4,64);
x = linspace( 0,4,64);
[SS,XX,FF,FAM] = compute_fa_2D(data,'i',s,x,T);

% FF = FF - FAM;

figure
contourf(SS,XX,(FF),20)
% pca = pcolor(SS,XX,(FF)); set(pca,'EdgeColor','None');
% set(gca,'ColorScale','log'); %shading interp;
% clim([1e-4 1]);
hold on
% plot(s,(s.^2/(1+data.fort_00.GEOMETRY.eps)),'--w');
% plot(s,eps_eff*(s.^2),'--k');
colormap(bluewhitered)