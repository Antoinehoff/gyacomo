DATA.CO      = h5readatt(filename,'/data/input/coll','CO');
DATA.K_N     = h5readatt(filename,'/data/input/ions','k_N');
DATA.K_T     = h5readatt(filename,'/data/input/ions','k_T');
DATA.Q0      = h5readatt(filename,'/data/input/geometry','q0');
DATA.EPS     = h5readatt(filename,'/data/input/geometry','eps');
DATA.SHEAR   = h5readatt(filename,'/data/input/geometry','shear');
DATA.GEOM    = h5readatt(filename,'/data/input/geometry','geometry');
% DATA.KAPPA   = h5readatt(filename,'/data/input/geometry','kappa');
% DATA.DELTA   = h5readatt(filename,'/data/input/geometry','delta');

DATA.DT_SIM  = h5readatt(filename,'/data/input/basic','dt');
DATA.PMAX    = h5readatt(filename,'/data/input/grid','pmax');
DATA.JMAX    = h5readatt(filename,'/data/input/grid','jmax');
DATA.Nx      = h5readatt(filename,'/data/input/grid','Nx');
DATA.Ny      = h5readatt(filename,'/data/input/grid','Ny');
DATA.L       = h5readatt(filename,'/data/input/grid','Lx');
DATA.CLOS    = h5readatt(filename,'/data/input/model','CLOS');
DATA.NL_CLOS = h5readatt(filename,'/data/input/model','NL_CLOS');
DATA.Na      = h5readatt(filename,'/data/input/model','Na');
DATA.NU      = h5readatt(filename,'/data/input/model','nu');
DATA.MUp     = h5readatt(filename,'/data/input/model','mu_p');
DATA.MUj     = h5readatt(filename,'/data/input/model','mu_j');
DATA.MUx     = h5readatt(filename,'/data/input/model','mu_x');
DATA.MUy     = h5readatt(filename,'/data/input/model','mu_y');
DATA.MUz     = h5readatt(filename,'/data/input/model','mu_z');
DATA.LINEARITY = h5readatt(filename,'/data/input/model','LINEARITY');
DATA.BETA    = h5readatt(filename,'/data/input/model','beta');
DATA.TAU_E   = h5readatt(filename,'/data/input/model','tau_e');
DATA.HYP_V   = h5readatt(filename,'/data/input/model','HYP_V');
DATA.K_cB    = h5readatt(filename,'/data/input/model','k_cB');
DATA.K_gB    = h5readatt(filename,'/data/input/model','k_gB');

DATA.W_GAMMA   = h5readatt(filename,'/data/input/diag_par','write_gamma') == 'y';
DATA.W_PHI     = h5readatt(filename,'/data/input/diag_par','write_phi')   == 'y';
DATA.W_NA00    = h5readatt(filename,'/data/input/diag_par','write_Na00')  == 'y';
DATA.W_NAPJ    = h5readatt(filename,'/data/input/diag_par','write_Napj')  == 'y';
DATA.W_SAPJ    = h5readatt(filename,'/data/input/diag_par','write_Sapj')  == 'y';

% Species dependent parameters
DATA.sigma = zeros(1,DATA.Na);
DATA.tau   = zeros(1,DATA.Na);
DATA.q     = zeros(1,DATA.Na);
DATA.K_N   = zeros(1,DATA.Na);
DATA.K_T   = zeros(1,DATA.Na);
spnames = {'ions','electrons'};
for ia=1:DATA.Na
    spdata = ['/data/input/',spnames{ia}];
    DATA.sigma(ia) = h5readatt(filename,spdata,'sigma');
    DATA.tau(ia)   = h5readatt(filename,spdata,'tau');
    DATA.q(ia)     = h5readatt(filename,spdata,'q');
    DATA.K_N(ia)   = h5readatt(filename,spdata,'k_N');
    DATA.K_T(ia)   = h5readatt(filename,spdata,'k_T');
end
DATA.spnames = spnames{1:DATA.Na};
DATA.CONAME = DATA.CO;

if    (DATA.CLOS == 0); DATA.CLOSNAME = 'Trunc.';
elseif(DATA.CLOS == 1); DATA.CLOSNAME = 'Clos. 1';
elseif(DATA.CLOS == 2); DATA.CLOSNAME = 'Clos. 2';
end

degngrad   = ['P_',num2str(DATA.PMAX),'_J_',num2str(DATA.JMAX)];

degngrad = [degngrad,'_Kni_%1.1f_nu_%0.0e_',...
        DATA.CONAME,'_CLOS_',num2str(DATA.CLOS),'_mu_%0.0e'];
degngrad   = sprintf(degngrad,[DATA.K_N,DATA.NU,DATA.MUx]);
% if ~DATA.LINEARITY; degngrad = ['lin_',degngrad]; end
resolution = [num2str(DATA.Nx),'x',num2str(DATA.Ny),'_'];
gridname   = ['L_',num2str(DATA.L),'_'];
DATA.PARAMS = [resolution,gridname,degngrad];