%% Auxiliary script to save figure using a dir (FIGDIR), name (FIGNAME) 
%  and parameters
FIGDIR = ['../results/', BASIC.SIMID,'/'];
if ~exist(FIGDIR, 'dir')
   mkdir(FIGDIR)
end

FIGNAME = [FIGDIR, FIGNAME,FMT];
saveas(fig,FIGNAME);
disp(['Figure saved @ : ',FIGNAME])