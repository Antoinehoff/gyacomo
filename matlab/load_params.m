CO      = h5readatt(filename,'/data/input','CO');
ETAB    = h5readatt(filename,'/data/input','eta_B');
ETAN    = h5readatt(filename,'/data/input','eta_n');
ETAT    = h5readatt(filename,'/data/input','eta_T');
PMAXI   = h5readatt(filename,'/data/input','pmaxi');
JMAXI   = h5readatt(filename,'/data/input','jmaxi');
PMAXE   = h5readatt(filename,'/data/input','pmaxe');
JMAXE   = h5readatt(filename,'/data/input','jmaxe');
NON_LIN = h5readatt(filename,'/data/input','NON_LIN');
NU      = h5readatt(filename,'/data/input','nu')/0.532;
NR      = h5readatt(filename,'/data/input','nr');
NZ      = h5readatt(filename,'/data/input','nz');
L       = h5readatt(filename,'/data/input','Lr');
if NON_LIN == 'y'
    NON_LIN = 1;
else
    NON_LIN = 0;
end

if    (CO == 0); CONAME = 'LB';
elseif(CO == -1); CONAME = 'FC';
elseif(CO == -2); CONAME = 'DG';
elseif(CO == -3); CONAME = 'DGGK';
end

if    (CLOS == 0); CLOSNAME = 'Trunc.';
elseif(CLOS == 1); CLOSNAME = 'Clos. 1';
elseif(CLOS == 2); CLOSNAME = 'Clos. 2';
end
if (PMAXE == PMAXI) && (JMAXE == JMAXI)
    degngrad   = ['P_',num2str(PMAXE),'_J_',num2str(JMAXE)];
else
    degngrad   = ['Pe_',num2str(PMAXE),'_Je_',num2str(JMAXE),...
        '_Pi_',num2str(PMAXI),'_Ji_',num2str(JMAXI)];
end
degngrad = [degngrad,'_eta_',num2str(ETAB/ETAN),'_nu_%0.0e_',...
        CONAME,'_CLOS_',num2str(CLOS),'_mu_%0.0e'];
degngrad   = sprintf(degngrad,[NU,MU]);
if ~NON_LIN; degngrad = ['lin_',degngrad]; end
resolution = [num2str(NR),'x',num2str(NZ/2),'_'];
gridname   = ['L_',num2str(L),'_'];
PARAMS = [resolution,gridname,degngrad];
% BASIC.RESDIR = [SIMDIR,PARAMS,'/'];
