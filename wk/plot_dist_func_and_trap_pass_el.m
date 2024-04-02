options.SHOW_FLUXSURF = 0;
options.SHOW_METRICS  = 0;
options.SHOW_CURVOP   = 0;
[fig, geo_arrays] = plot_metric(data,options);

eps_eff = 1/max((geo_arrays.hatB));
T = [10];
s = linspace(-3,3,64);
x = linspace(0,8,48);
[SS,XX,FF] = compute_fa_2D(data,'i',s,x,T);

figure
% contour(SS,XX,log10(FF))
pca = pcolor(SS,XX,(FF)); set(pca,'EdgeColor','None');
set(gca,'ColorScale','log'); %shading interp;
% clim([1e-4 1]);
hold on
plot(s,(s.^2/(1+data.fort_00.GEOMETRY.eps)),'--w');
plot(s,eps_eff*(s.^2),'--k');