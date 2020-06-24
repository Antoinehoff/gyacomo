%% Auxiliary script to save figure using a dir (FIGDIR), name (FIGNAME) 
%  and parameters
FIGDIR = ['../results/', SIMID,'/'];
if ~exist(FIGDIR, 'dir')
   mkdir(FIGDIR)
end
FIGNAME = [FIGNAME,'_Pe_',num2str(GRID.pmaxe),'_Je_',num2str(GRID.jmaxe),...
    '_Pi_',num2str(GRID.pmaxi),'_Ji_',num2str(GRID.jmaxi),...
    '_etan_',num2str(MODEL.eta_n),'_etaB_',num2str(MODEL.eta_B),'_nu_',num2str(MODEL.nu)];
FIGNAME = [FIGDIR, FIGNAME,'.fig'];
savefig(fig,FIGNAME);
disp(['Figure saved @ : ',FIGNAME])