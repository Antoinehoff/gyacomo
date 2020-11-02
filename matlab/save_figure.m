%% Auxiliary script to save figure using a dir (FIGDIR), name (FIGNAME) 
%  and parameters
if ~exist([BASIC.RESDIR,'/fig'], 'dir')
   mkdir([BASIC.RESDIR,'/fig'])
end
saveas(fig,[BASIC.RESDIR,'/fig/', FIGNAME,'.fig']);
saveas(fig,[BASIC.RESDIR, FIGNAME,'.png']);
disp(['Figure saved @ : ',[BASIC.RESDIR, FIGNAME,'.png']])