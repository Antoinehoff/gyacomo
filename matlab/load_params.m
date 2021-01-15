CONAME = h5readatt(filename,'/data/input','CO');
if CONAME == -1
    CONAME = 'FC';
elseif CONAME == 0
    CONAME = 'LB';
elseif CONAME == -2
    CONAME = 'DG';
elseif CONAME == -3
    CONAME = 'DGGK';
end
ETAB = h5readatt(filename,'/data/input','eta_B');
ETAN = h5readatt(filename,'/data/input','eta_n');
ETAT = h5readatt(filename,'/data/input','eta_T');
PMAXI = h5readatt(filename,'/data/input','pmaxi');
JMAXI = h5readatt(filename,'/data/input','jmaxi');
PMAXE = h5readatt(filename,'/data/input','pmaxe');
JMAXE = h5readatt(filename,'/data/input','jmaxe');
NON_LIN = h5readatt(filename,'/data/input','NON_LIN');
NU = h5readatt(filename,'/data/input','nu')/0.532;
if NON_LIN == 'y'
    NON_LIN = 1;
else
    NON_LIN = 0;
end

degngrad   = ['Pe_',num2str(PMAXE),'_Je_',num2str(JMAXE),...
    '_Pi_',num2str(PMAXI),'_Ji_',num2str(JMAXI),...
    '_nB_',num2str(ETAB),'_nN_',num2str(ETAN),'_nu_%0.0e_',...
    CONAME,'_mu_%0.0e'];
degngrad   = sprintf(degngrad,[NU,MU]);
if ~NON_LIN; degngrad = ['lin_',degngrad]; end
resolution = [num2str(GRID.Nr),'x',num2str(GRID.Nz/2),'_'];
gridname   = ['L_',num2str(L),'_'];
PARAMS = [resolution,gridname,degngrad];
BASIC.RESDIR = [SIMDIR,PARAMS,'/'];
