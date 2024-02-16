function [DATA] = load_params(DATA,filename)
    DATA.inputs.CO      = h5readatt(filename,'/data/input/coll','CO');

    DATA.inputs.Q0      = h5readatt(filename,'/data/input/geometry','q0');
    DATA.inputs.EPS     = h5readatt(filename,'/data/input/geometry','eps');
    DATA.inputs.SHEAR   = h5readatt(filename,'/data/input/geometry','shear');
    DATA.inputs.GEOM    = h5readatt(filename,'/data/input/geometry','geometry');
    % DATA.KAPPA   = h5readatt(filename,'/data/input/geometry','kappa');
    % DATA.DELTA   = h5readatt(filename,'/data/input/geometry','delta');

    DATA.inputs.DT_SIM  = h5readatt(filename,'/data/input/basic','dt');
    DATA.inputs.PMAX    = h5readatt(filename,'/data/input/grid','pmax');
    DATA.inputs.JMAX    = h5readatt(filename,'/data/input/grid','jmax');
    DATA.inputs.Nx      = h5readatt(filename,'/data/input/grid','Nx');
    DATA.inputs.Ny      = h5readatt(filename,'/data/input/grid','Ny');
    DATA.inputs.L       = h5readatt(filename,'/data/input/grid','Lx');
    try
    DATA.inputs.CLOS    = h5readatt(filename,'/data/input/model','CLOS');
    DATA.inputs.NL_CLOS = h5readatt(filename,'/data/input/model','NL_CLOS');
    catch
        try
            DATA.inputs.ha_cl   = h5readatt(filename,'/data/input/closure','hierarchy_closure');
            DATA.inputs.CLOS    = h5readatt(filename,'/data/input/closure','dmax');   
            DATA.inputs.nl_cl   = h5readatt(filename,'/data/input/closure','nonlinear_closure');   
            DATA.inputs.NL_CLOS = h5readatt(filename,'/data/input/closure','nmax');   
        catch
            DATA.inputs.CLOS = 99;
            DATA.inputs.NL_CLOS = 99;
        end
    end
    DATA.inputs.Na      = h5readatt(filename,'/data/input/model','Na');
    DATA.inputs.NU      = h5readatt(filename,'/data/input/model','nu');
    DATA.inputs.MUp     = h5readatt(filename,'/data/input/model','mu_p');
    DATA.inputs.MUj     = h5readatt(filename,'/data/input/model','mu_j');
    DATA.inputs.MUx     = h5readatt(filename,'/data/input/model','mu_x');
    DATA.inputs.MUy     = h5readatt(filename,'/data/input/model','mu_y');
    DATA.inputs.MUz     = h5readatt(filename,'/data/input/model','mu_z');
    DATA.inputs.LINEARITY = h5readatt(filename,'/data/input/model','LINEARITY');
    DATA.inputs.BETA    = h5readatt(filename,'/data/input/model','beta');
    try
        DATA.inputs.TAU_I   = h5readatt(filename,'/data/input/model','tau_i');
    catch
        DATA.inputs.TAU_I   = 1.0;
    end
    DATA.inputs.HYP_V   = h5readatt(filename,'/data/input/model','HYP_V');
    DATA.inputs.K_cB    = h5readatt(filename,'/data/input/model','k_cB');
    DATA.inputs.K_gB    = h5readatt(filename,'/data/input/model','k_gB');
    DATA.inputs.ADIAB_E = h5readatt(filename,'/data/input/model','ADIAB_E') == 'y';
    try
        DATA.inputs.ADIAB_I = h5readatt(filename,'/data/input/model','ADIAB_I') == 'y';
    catch 
        DATA.inputs.ADIAB_I = 0;
    end
    try
        DATA.inputs.W_GAMMA   = h5readatt(filename,'/data/input/diagnostics','write_gamma') == 'y';
        DATA.inputs.W_PHI     = h5readatt(filename,'/data/input/diagnostics','write_phi')   == 'y';
        DATA.inputs.W_NA00    = h5readatt(filename,'/data/input/diagnostics','write_Na00')  == 'y';
        DATA.inputs.W_NAPJ    = h5readatt(filename,'/data/input/diagnostics','write_Napj')  == 'y';
        DATA.inputs.W_SAPJ    = h5readatt(filename,'/data/input/diagnostics','write_Sapj')  == 'y';
    catch
        DATA.inputs.W_GAMMA   = h5readatt(filename,'/data/input/diag_par','write_gamma') == 'y';
        DATA.inputs.W_PHI     = h5readatt(filename,'/data/input/diag_par','write_phi')   == 'y';
        DATA.inputs.W_NA00    = h5readatt(filename,'/data/input/diag_par','write_Na00')  == 'y';
        DATA.inputs.W_NAPJ    = h5readatt(filename,'/data/input/diag_par','write_Napj')  == 'y';
        DATA.inputs.W_SAPJ    = h5readatt(filename,'/data/input/diag_par','write_Sapj')  == 'y';
    end
    % Species dependent parameters
    DATA.inputs.sigma = zeros(1,DATA.inputs.Na);
    DATA.inputs.tau   = zeros(1,DATA.inputs.Na);
    DATA.inputs.q     = zeros(1,DATA.inputs.Na);
    DATA.inputs.K_N   = zeros(1,DATA.inputs.Na);
    DATA.inputs.K_T   = zeros(1,DATA.inputs.Na);
    spnames = {'ions','electrons'};
    if(DATA.inputs.ADIAB_E) 
        spnames = {spnames{1}};
    end
    if(DATA.inputs.ADIAB_I) 
        spnames = {spnames{2}};
    end
        
    for ia=1:DATA.inputs.Na
        spdata = ['/data/input/',spnames{ia}];
        DATA.inputs.sigma(ia) = h5readatt(filename,spdata,'sigma');
        DATA.inputs.tau(ia)   = h5readatt(filename,spdata,'tau');
        DATA.inputs.q(ia)     = h5readatt(filename,spdata,'q');
        DATA.inputs.K_N(ia)   = h5readatt(filename,spdata,'k_N');
        DATA.inputs.K_T(ia)   = h5readatt(filename,spdata,'k_T');
    end
    DATA.inputs.spnames = spnames{1:DATA.inputs.Na};
    DATA.inputs.CONAME = DATA.inputs.CO;

    if    (DATA.inputs.CLOS == 0); DATA.CLOSNAME = 'Trunc.';
    elseif(DATA.inputs.CLOS == 1); DATA.CLOSNAME = 'Clos. 1';
    elseif(DATA.inputs.CLOS == 2); DATA.CLOSNAME = 'Clos. 2';
    end

    degngrad   = ['P_',num2str(DATA.inputs.PMAX),'_J_',num2str(DATA.inputs.JMAX)];

    degngrad = [degngrad,'_Kni_%1.1f_nu_%0.0e_',...
            DATA.inputs.CONAME,'_CLOS_',num2str(DATA.inputs.CLOS),'_mu_%0.0e'];
    degngrad   = sprintf(degngrad,[DATA.inputs.K_N,DATA.inputs.NU,DATA.inputs.MUx]);
    % if ~DATA.LINEARITY; degngrad = ['lin_',degngrad]; end
    resolution = [num2str(DATA.inputs.Nx),'x',num2str(DATA.inputs.Ny),'_'];
    gridname   = ['L_',num2str(DATA.inputs.L),'_'];
    DATA.params_string = [resolution,gridname,degngrad];
end