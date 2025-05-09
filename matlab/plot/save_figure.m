%% Auxiliary script to save figure using a dir (FIGDIR), name (FIGNAME) 
%  and parameters
function save_figure(data,FIGURE,EXT)
if ~exist([data.FIGDIR,'/fig'], 'dir')
   mkdir([data.FIGDIR,'/fig'])
end
saveas(FIGURE.fig,[data.FIGDIR,'/fig/', FIGURE.FIGNAME,'.fig']);
saveas(FIGURE.fig,[data.FIGDIR, FIGURE.FIGNAME,EXT]);
disp(['Figure saved @ : ',[data.FIGDIR, FIGURE.FIGNAME,EXT]])
end