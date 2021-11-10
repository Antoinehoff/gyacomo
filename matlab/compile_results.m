function [DATA] = compile_results(DIRECTORY,JOBNUMMIN,JOBNUMMAX)
DATA = {};
CONTINUE = 1;
JOBNUM   = JOBNUMMIN; JOBFOUND = 0;
DATA.TJOB_SE  = []; % Start and end times of jobs
DATA.NU_EVOL  = []; % evolution of parameter nu between jobs
DATA.CO_EVOL  = []; % evolution of CO
DATA.MU_EVOL  = []; % evolution of parameter mu between jobs
DATA.K_N_EVOL = []; %
DATA.L_EVOL   = []; % 
DATA.DT_EVOL  = []; %
% FIELDS
Nipj_    = []; Nepj_    = [];
Ni00_    = []; Ne00_    = [];
HFLUX_   = [];
GGAMMA_  = [];
PGAMMA_  = [];
PHI_     = [];
DENS_E_  = [];
DENS_I_  = [];
TEMP_E_  = [];
TEMP_I_  = [];
Ts0D_    = [];
Ts3D_    = [];
Ts5D_    = [];
Sipj_    = []; Sepj_    = [];
Pe_old   = 1e9; Pi_old = Pe_old; Je_old = Pe_old; Ji_old = Pe_old;
Pi_max=0; Pe_max=0; Ji_max=0; Je_max=0;
while(CONTINUE)
    filename = sprintf([DIRECTORY,'outputs_%.2d.h5'],JOBNUM);
    if (exist(filename, 'file') == 2 && JOBNUM <= JOBNUMMAX)
        % Load results of simulation #JOBNUM
%         load_results
        %% load results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['Loading ',filename])
        % Loading from output file
        CPUTIME   = h5readatt(filename,'/data/input','cpu_time');
        DT_SIM    = h5readatt(filename,'/data/input','dt');
        [Pe, Je, Pi, Ji, kx, ky, z] = load_grid_data(filename);
        W_GAMMA   = strcmp(h5readatt(filename,'/data/input','write_gamma'),'y');
        W_HF      = 0;%strcmp(h5readatt(filename,'/data/input','write_hf'   ),'y');
        W_PHI     = strcmp(h5readatt(filename,'/data/input','write_phi'  ),'y');
        W_NA00    = strcmp(h5readatt(filename,'/data/input','write_Na00' ),'y');
        W_NAPJ    = strcmp(h5readatt(filename,'/data/input','write_Napj' ),'y');
        W_SAPJ    = strcmp(h5readatt(filename,'/data/input','write_Sapj' ),'y');
        W_DENS    = strcmp(h5readatt(filename,'/data/input','write_dens' ),'y');
        W_TEMP    = strcmp(h5readatt(filename,'/data/input','write_temp' ),'y');
        KIN_E     = strcmp(h5readatt(filename,'/data/input',     'KIN_E' ),'y');
%         KIN_E     = 1;
        
        % Check polynomials degrees
        Pe_new= numel(Pe); Je_new= numel(Je);
        Pi_new= numel(Pi); Ji_new= numel(Ji);
        if(Pe_max < Pe_new); Pe_max = Pe_new; end;
        if(Je_max < Je_new); Je_max = Je_new; end;
        if(Pi_max < Pi_new); Pi_max = Pi_new; end;
        if(Ji_max < Ji_new); Ji_max = Ji_new; end;
        % If a degree is larger than previous job, put them in a larger array
        if (sum([Pe_new, Je_new, Pi_new, Ji_new]>[Pe_old, Je_old, Pi_old, Ji_old]) >= 1)
            if W_NAPJ
                tmp = Nipj_; sz = size(tmp);
                Nipj_ = zeros(cat(1,[Pi_new,Ji_new]',sz(3:end)')');
                Nipj_(1:Pi_old,1:Ji_old,:,:,:) = tmp;
                tmp = Nepj_; sz = size(tmp);
                Nepj_ = zeros(cat(1,[Pe_new,Je_new]',sz(3:end)')');
                Nepj_(1:Pe_old,1:Je_old,:,:,:) = tmp;
            end
            if W_SAPJ
                tmp = Sipj_; sz = size(tmp);
                Sipj_ = zeros(cat(1,[Pi_new,Ji_new]',sz(3:end)')');
                Sipj_(1:Pi_old,1:Ji_old,:,:,:) = tmp;
                tmp = Sepj_; sz = size(tmp);
                Sepj_ = zeros(cat(1,[Pe_new,Je_new]',sz(3:end)')');
                Sepj_(1:Pe_old,1:Je_old,:,:,:) = tmp;
            end
        % If a degree is smaller than previous job, put zero to add. deg.
        elseif (sum([Pe_new, Je_new, Pi_new, Ji_new]<[Pe_old, Je_old, Pi_old, Ji_old]) >= 1 && Pe_old ~= 1e9)
            if W_NAPJ
                tmp = Nipj; sz = size(tmp);
                Nipj = zeros(cat(1,[Pi_max,Ji_max]',sz(3:end)')');
                Nipj(1:Pi_new,1:Ji_new,:,:,:) = tmp;
                tmp = Nepj; sz = size(tmp);
                Nepj = zeros(cat(1,[Pe_max,Je_max]',sz(3:end)')');
                Nepj(1:Pe_new,1:Je_new,:,:,:) = tmp;
            end
            if W_SAPJ
                tmp = Sipj; sz = size(tmp);
                Sipj = zeros(cat(1,[Pi_max,Ji_max]',sz(3:end)')');
                Sipj(1:Pi_new,1:Ji_new,:,:,:) = tmp;
                tmp = Sepj; sz = size(tmp);
                Sepj = zeros(cat(1,[Pe_max,Je_max]',sz(3:end)')');
                Sepj(1:Pe_new,1:Je_new,:,:,:) = tmp;
            end
        end


        if W_GAMMA
            [ GGAMMA_RI, Ts0D, ~] = load_0D_data(filename, 'gflux_ri');
            PGAMMA_RI            = load_0D_data(filename, 'pflux_ri');
            GGAMMA_ = cat(1,GGAMMA_,GGAMMA_RI); clear GGAMMA_RI
            PGAMMA_ = cat(1,PGAMMA_,PGAMMA_RI); clear PGAMMA_RI
        end

        if W_HF
            [ HFLUX_X, Ts0D, ~] = load_0D_data(filename, 'hflux_x');
            HFLUX_ = cat(1,HFLUX_,HFLUX_X); clear HFLUX_X
        end

        if W_PHI
            [ PHI, Ts3D, ~] = load_3D_data(filename, 'phi');
            PHI_  = cat(4,PHI_,PHI); clear PHI
        end
        if W_NA00
            if KIN_E
             Ne00  = load_3D_data(filename, 'Ne00');
             Ne00_ = cat(4,Ne00_,Ne00); clear Ne00
            end
            [Ni00, Ts3D, ~] = load_3D_data(filename, 'Ni00');
             Ni00_ = cat(4,Ni00_,Ni00); clear Ni00
        end
        if W_DENS
            if KIN_E
        [DENS_E, Ts3D, ~] = load_3D_data(filename, 'dens_e');
            DENS_E_ = cat(4,DENS_E_,DENS_E); clear DENS_E
            end
        [DENS_I, Ts3D, ~] = load_3D_data(filename, 'dens_i');
            DENS_I_ = cat(4,DENS_I_,DENS_I); clear DENS_I
        end
        if W_TEMP
            if KIN_E 
        [TEMP_E, Ts3D, ~] = load_3D_data(filename, 'temp_e');
                TEMP_E_ = cat(4,TEMP_E_,TEMP_E); clear TEMP_E
            end
        [TEMP_I, Ts3D, ~] = load_3D_data(filename, 'temp_i');
            TEMP_I_ = cat(4,TEMP_I_,TEMP_I); clear TEMP_I
        end

        if W_NAPJ
        [Nipj, Ts5D, ~] = load_5D_data(filename, 'moments_i');
            Nipj_ = cat(6,Nipj_,Nipj); clear Nipj
            if KIN_E
                Nepj  = load_5D_data(filename, 'moments_e');
                Nepj_ = cat(6,Nepj_,Nepj); clear Nepj
            end
        end
        if W_SAPJ
     	  Sipj_ = cat(6,Sipj_,Sipj);
          if KIN_E
           Sepj_ = cat(6,Sepj_,Sepj);
          end
        end
        Ts0D_   = cat(1,Ts0D_,Ts0D);
        Ts3D_   = cat(1,Ts3D_,Ts3D);
        Ts5D_   = cat(1,Ts5D_,Ts5D);

        % Evolution of simulation parameters
        load_params
        DATA.TJOB_SE   = [DATA.TJOB_SE  Ts0D(1) Ts0D(end)];
        DATA.NU_EVOL   = [DATA.NU_EVOL  DATA.NU     DATA.NU];
        DATA.CO_EVOL   = [DATA.CO_EVOL  DATA.CO     DATA.CO];
        DATA.MU_EVOL   = [DATA.MU_EVOL  DATA.MU     DATA.MU];
        DATA.K_N_EVOL  = [DATA.K_N_EVOL DATA.K_N    DATA.K_N];
        DATA.L_EVOL    = [DATA.L_EVOL   DATA.L      DATA.L];
        DATA.DT_EVOL   = [DATA.DT_EVOL  DATA.DT_SIM DATA.DT_SIM];

        JOBFOUND = JOBFOUND + 1;
        LASTJOB  = JOBNUM;
        Pe_old = Pe_new; Je_old = Je_new;
        Pi_old = Pi_new; Ji_old = Ji_new;
    elseif (JOBNUM > JOBNUMMAX)
        CONTINUE = 0;
        disp(['found ',num2str(JOBFOUND),' results']);
    end
    JOBNUM   = JOBNUM + 1;
end

%% Build grids
Nkx = numel(kx); Nky = numel(ky);
[KY,KX] = meshgrid(ky,kx);
Lkx = max(kx)-min(kx); Lky = max(ky)-min(ky);
dkx = Lkx/(Nkx-1); dky = Lky/(Nky-1);
KPERP2 = KY.^2+KX.^2;
[~,ikx0] = min(abs(kx)); [~,iky0] = min(abs(ky));

Nx = 2*(Nkx-1);  Ny = Nky;      Nz = numel(z);
Lx = 2*pi/dkx;   Ly = 2*pi/dky;
dx = Lx/Nx;      dy = Ly/Ny; dz = 2*pi/Nz;
x = dx*(-Nx/2:(Nx/2-1)); Lx = max(x)-min(x);
y = dy*(-Ny/2:(Ny/2-1)); Ly = max(y)-min(y);
%% Add everything in output structure
% Fields
DATA.GGAMMA_RI = GGAMMA_; DATA.PGAMMA_RI = PGAMMA_; DATA.HFLUX_X = HFLUX_;
DATA.Nipj = Nipj_; DATA.Ni00 = Ni00_; DATA.DENS_I = DENS_I_; DATA.TEMP_I = TEMP_I_;
if(KIN_E)
DATA.Nepj = Nepj_; DATA.Ne00 = Ne00_; DATA.DENS_E = DENS_E_; DATA.TEMP_E = TEMP_E_;
end
DATA.Ts5D = Ts5D_; DATA.Ts3D = Ts3D_; DATA.Ts0D = Ts0D_;
DATA.PHI  = PHI_; 
% grids
DATA.Pe = Pe; DATA.Pi = Pi; 
DATA.kx = kx; DATA.ky = ky; DATA.z = z;
DATA.x  = x;  DATA.y  = y;
DATA.ikx0 = ikx0; DATA.iky0 = iky0;
DATA.Nx = Nx; DATA.Ny = Ny; DATA.Nz = Nz; DATA.Nkx = Nkx; DATA.Nky = Nky; 
DATA.Pmaxe = numel(Pe); DATA.Pmaxi = numel(Pi); DATA.Jmaxe = numel(Je); DATA.Jmaxi = numel(Ji);
DATA.dir      = DIRECTORY;
DATA.localdir = ['..',DIRECTORY(20:end)];

JOBNUM = LASTJOB;

filename = sprintf([DIRECTORY,'outputs_%.2d.h5'],JOBNUM);
end