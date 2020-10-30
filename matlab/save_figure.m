%% Auxiliary script to save figure using a dir (FIGDIR), name (FIGNAME) 
%  and parameters
FIGNAME = [BASIC.RESDIR, FIGNAME,FMT];
saveas(fig,FIGNAME);
disp(['Figure saved @ : ',FIGNAME])