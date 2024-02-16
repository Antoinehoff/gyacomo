options.NAME      = ['N_i^{00}'];
options.PLAN      = 'xyz';
options.TIME      = [30];
TOPLOT = process_field( data, options );

[X,Y,Z] = meshgrid(data.grids.x,data.grids.y,data.grids.z);
X = smooth3(X,'gaussian',5);
Y = smooth3(Y,'gaussian',5);
Z = smooth3(Z,'gaussian',5);
TOPLOT.FIELD = smooth3(TOPLOT.FIELD,'gaussian',5);
s = isosurface(X,Y,Z,TOPLOT.FIELD,0);
s = reducepatch(s,0.1);
%% figure

figure

p=patch(s);
set(p,'EdgeColor','none');
set(p,'FaceColor','r');
set(p,'FaceAlpha',0.4);
