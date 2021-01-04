CONAME = h5readatt(filename,'/data/input','CO');
if CONAME == -1
    CONAME = 'FC';
elseif CONAME == 0
    CONAME = 'LB';
elseif CONAME == 1
    CONAME = 'DG';
end
ETAB = h5readatt(filename,'/data/input','eta_B');
ETAN = h5readatt(filename,'/data/input','eta_n');
ETAT = h5readatt(filename,'/data/input','eta_T');

NON_LIN = h5readatt(filename,'/data/input','NON_LIN');
if NON_LIN == 'y'
    NON_LIN = 1;
else
    NON_LIN = 0;
end
