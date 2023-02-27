DATA.CO      = h5readatt(filename,'/data/input','CO');
try
    DATA.ETA_N   = h5readatt(filename,'/data/input','ETA_N');
    DATA.ETA_T   = h5readatt(filename,'/data/input','ETA_T');
catch
    DATA.ETA_N = 1.0;
    DATA.ETA_T = 1.0;
end
try
    DATA.K_N     = h5readatt(filename,'/data/input','K_n');
catch
    try
        DATA.K_N     = h5readatt(filename,'/data/input','k_N');
    catch
        DATA.K_N     = h5readatt(filename,'/data/input','k_Ni');
    end
end
try
    DATA.K_T     = h5readatt(filename,'/data/input','K_T');
catch
    try
        DATA.K_T     = h5readatt(filename,'/data/input','k_T');
    catch
        DATA.K_T     = h5readatt(filename,'/data/input','k_Ti');
    end
end
DATA.sigma_e = h5readatt(filename,'/data/input','sigma_e');
DATA.sigma_i = h5readatt(filename,'/data/input','sigma_i');
DATA.tau_e   = h5readatt(filename,'/data/input','tau_e');
DATA.tau_i   = h5readatt(filename,'/data/input','tau_i');
DATA.q_e     = h5readatt(filename,'/data/input','q_e');
DATA.q_i     = h5readatt(filename,'/data/input','q_i');
DATA.Q0      = h5readatt(filename,'/data/input','q0');
DATA.SHEAR   = h5readatt(filename,'/data/input','shear');
DATA.EPS     = h5readatt(filename,'/data/input','eps');
DATA.PMAXI   = h5readatt(filename,'/data/input','pmaxi');
DATA.JMAXI   = h5readatt(filename,'/data/input','jmaxi');
DATA.PMAXE   = h5readatt(filename,'/data/input','pmaxe');
DATA.JMAXE   = h5readatt(filename,'/data/input','jmaxe');
% DATA.LINEARITY = h5readatt(filename,'/data/input','NON_LIN');
% DATA.LINEARITY = h5readatt(filename,'/data/input','LINEARITY');
DATA.NU      = h5readatt(filename,'/data/input','nu');
DATA.Nx      = h5readatt(filename,'/data/input','Nx');
DATA.Ny      = h5readatt(filename,'/data/input','Ny');
DATA.L       = h5readatt(filename,'/data/input','Lx');
DATA.CLOS    = h5readatt(filename,'/data/input','CLOS');
DATA.DT_SIM  = h5readatt(filename,'/data/input','dt');
DATA.MU      = h5readatt(filename,'/data/input','mu');
DATA.MUx     = h5readatt(filename,'/data/input','mu_x');
DATA.MUy     = h5readatt(filename,'/data/input','mu_y');
DATA.MUz     = h5readatt(filename,'/data/input','mu_z');
try
    DATA.BETA    = h5readatt(filename,'/data/input','beta');
catch
    DATA.BETA  = 0;
end
DATA.W_GAMMA   = h5readatt(filename,'/data/input','write_gamma') == 'y';
DATA.W_PHI     = h5readatt(filename,'/data/input','write_phi')   == 'y';
DATA.W_NA00    = h5readatt(filename,'/data/input','write_Na00')  == 'y';
DATA.W_NAPJ    = h5readatt(filename,'/data/input','write_Napj')  == 'y';
DATA.W_SAPJ    = h5readatt(filename,'/data/input','write_Sapj')  == 'y';

% if DATA.LINEARITY == 'y'
%     DATA.LINEARITY = 1;
% else
%     DATA.LINEARITY = 0;
% end

DATA.CONAME = DATA.CO;

if    (DATA.CLOS == 0); DATA.CLOSNAME = 'Trunc.';
elseif(DATA.CLOS == 1); DATA.CLOSNAME = 'Clos. 1';
elseif(DATA.CLOS == 2); DATA.CLOSNAME = 'Clos. 2';
end
if (DATA.PMAXE == DATA.PMAXI) && (DATA.JMAXE == DATA.JMAXI)
    degngrad   = ['P_',num2str(DATA.PMAXE),'_J_',num2str(DATA.JMAXE)];
else
    degngrad   = ['Pe_',num2str(DATA.PMAXE),'_Je_',num2str(DATA.JMAXE),...
        '_Pi_',num2str(DATA.PMAXI),'_Ji_',num2str(DATA.JMAXI)];
end
degngrad = [degngrad,'_Kni_%1.1f_nu_%0.0e_',...
        DATA.CONAME,'_CLOS_',num2str(DATA.CLOS),'_mu_%0.0e'];
degngrad   = sprintf(degngrad,[DATA.K_N,DATA.NU,DATA.MU]);
% if ~DATA.LINEARITY; degngrad = ['lin_',degngrad]; end
resolution = [num2str(DATA.Nx),'x',num2str(DATA.Ny),'_'];
gridname   = ['L_',num2str(DATA.L),'_'];
DATA.PARAMS = [resolution,gridname,degngrad];