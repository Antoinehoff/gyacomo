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
