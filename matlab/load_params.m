CO      = h5readatt(filename,'/data/input','CO');
ETAB    = h5readatt(filename,'/data/input','eta_B');
ETAN    = h5readatt(filename,'/data/input','eta_n');
ETAT    = h5readatt(filename,'/data/input','eta_T');
PMAXI   = h5readatt(filename,'/data/input','pmaxi');
JMAXI   = h5readatt(filename,'/data/input','jmaxi');
PMAXE   = h5readatt(filename,'/data/input','pmaxe');
JMAXE   = h5readatt(filename,'/data/input','jmaxe');
NON_LIN = h5readatt(filename,'/data/input','NON_LIN');
NU      = h5readatt(filename,'/data/input','nu');
Nx      = h5readatt(filename,'/data/input','Nx');
Ny      = h5readatt(filename,'/data/input','Ny');
L       = h5readatt(filename,'/data/input','Lx');
CLOS    = h5readatt(filename,'/data/input','CLOS');
DT_SIM  = h5readatt(filename,'/data/input','dt');
MU      = h5readatt(filename,'/data/input','mu');
% MU      = str2num(filename(end-18:end-14)); %bad...
W_GAMMA   = h5readatt(filename,'/data/input','write_gamma') == 'y';
W_PHI     = h5readatt(filename,'/data/input','write_phi')   == 'y';
W_NA00    = h5readatt(filename,'/data/input','write_Na00')  == 'y';
W_NAPJ    = h5readatt(filename,'/data/input','write_Napj')  == 'y';
W_SAPJ    = h5readatt(filename,'/data/input','write_Sapj')  == 'y';

if NON_LIN == 'y'
    NON_LIN = 1;
else
    NON_LIN = 0;
end

if    (CO == -3); CONAME = 'PADK';
elseif(CO == -2); CONAME = 'SGDK';
elseif(CO == -1); CONAME = 'DGDK';
elseif(CO ==  0); CONAME = 'LB';
elseif(CO ==  1); CONAME = 'DGGK';
elseif(CO ==  2); CONAME = 'SGGK';
elseif(CO ==  3); CONAME = 'PAGK';
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
degngrad = [degngrad,'_eta_%1.1f_nu_%0.0e_',...
        CONAME,'_CLOS_',num2str(CLOS),'_mu_%0.0e'];
degngrad   = sprintf(degngrad,[ETAB/ETAN,NU,MU]);
if ~NON_LIN; degngrad = ['lin_',degngrad]; end
resolution = [num2str(Nx),'x',num2str(Ny/2),'_'];
gridname   = ['L_',num2str(L),'_'];
PARAMS = [resolution,gridname,degngrad];
% BASIC.RESDIR = [SIMDIR,PARAMS,'/'];
