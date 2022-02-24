%% Auxiliary script to save figure using a dir (FIGDIR), name (FIGNAME) 
%  and parameters
function save_figure(data,FIGURE)
if ~exist([data.localdir,'/fig'], 'dir')
   mkdir([data.localdir,'/fig'])
end
saveas(FIGURE.fig,[data.localdir,'/fig/', FIGURE.FIGNAME,'.fig']);
saveas(FIGURE.fig,[data.localdir, FIGURE.FIGNAME,'.png']);
disp(['Figure saved @ : ',[data.localdir, FIGURE.FIGNAME,'.png']])
end