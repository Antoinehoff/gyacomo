function [DATA] = compile_results_low_mem(DATA,DIRECTORY,JOBNUMMIN,JOBNUMMAX)
CONTINUE = 1;
JOBNUM   = JOBNUMMIN; JOBFOUND = 0;
DATA.CODENAME = 'GYAC';
DATA.folder = DIRECTORY;
DATA.params_evol.TJOB_SE  = []; % Start and end times of jobs
DATA.params_evol.NU_EVOL  = []; % evolution of parameter nu between jobs
DATA.params_evol.CO_EVOL  = []; % evolution of CO
DATA.params_evol.MUx_EVOL  = []; % evolution of parameter mu between jobs
DATA.params_evol.MUy_EVOL  = []; % evolution of parameter mu between jobs
DATA.params_evol.MUz_EVOL  = []; % evolution of parameter mu between jobs
DATA.params_evol.K_N_EVOL = []; %
DATA.params_evol.K_T_EVOL = []; %
DATA.params_evol.L_EVOL   = []; % 
DATA.params_evol.DT_EVOL  = []; %
% Low memoery cost data
HFLUX_   = [];
HFLUXE_   = [];
GGAMMA_ = [];
PGAMMA_ = [];
Ts0D_    = [];
DATA.outfilenames = {};
ii = 1;
while(CONTINUE)
    filename = sprintf([DIRECTORY,'outputs_%.2d.h5'],JOBNUM);
    % Check presence and jobnummax
    if (exist(filename, 'file') == 2 && JOBNUM <= JOBNUMMAX)
        DATA.outfilenames{ii} = filename;
        %test if it is corrupted or currently running
        try
            openable = ~isempty(h5read(filename,'/data/var3d/time'));
        catch
            openable = 0;
        end
        if openable
            %% load results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('Loading ID %.2d (%s)\n',JOBNUM,filename);
            % Loading from output file
            CPUTIME   = h5readatt(filename,'/data/input','cpu_time');
            DT_SIM    = h5readatt(filename,'/data/input/basic','dt');
            [Parray, Jarray, kx, ky, z] = load_grid_data(filename);
            W_GAMMA   = strcmp(h5readatt(filename,'/data/input/diagnostics','write_gamma'),'y');
            W_HF      = strcmp(h5readatt(filename,'/data/input/diagnostics','write_hf'   ),'y');
            KIN_E     = strcmp(h5readatt(filename,'/data/input/model',     'ADIAB_E' ),'n');
            BETA      = h5readatt(filename,'/data/input/model','beta');

            if W_GAMMA
                try
                    [ GGAMMA_R, Ts0D, ~] = load_0D_data(filename, 'gflux_x');
                    PGAMMA_R             = load_0D_data(filename, 'pflux_x');
                catch
                    % old version
                    [ GGAMMA_R, Ts0D, ~] = load_0D_data(filename, 'gflux_xi');
                    PGAMMA_R             = load_0D_data(filename, 'pflux_xi');
                end
                GGAMMA_ = cat(2,GGAMMA_,GGAMMA_R); clear GGAMMA_R
                PGAMMA_ = cat(2,PGAMMA_,PGAMMA_R); clear PGAMMA_R
            end

            if W_HF
                try
                    [ HFLUX_X, Ts0D, ~] = load_0D_data(filename, 'hflux_x');
                catch
                    % old version
                    [ HFLUX_X, Ts0D, ~] = load_0D_data(filename, 'hflux_xi');
                end
                HFLUX_ = cat(2,HFLUX_,HFLUX_X); clear HFLUX_X
            end

            Ts0D_   = cat(1,Ts0D_,Ts0D);

            % Evolution of simulation parameters
            DATA = load_params(DATA,filename);
            DATA.params_evol.TJOB_SE   = [DATA.params_evol.TJOB_SE   Ts0D(1)     Ts0D(end)];
            DATA.params_evol.NU_EVOL   = [DATA.params_evol.NU_EVOL   DATA.inputs.NU     DATA.inputs.NU];
            DATA.params_evol.CO_EVOL   = [DATA.params_evol.CO_EVOL   DATA.inputs.CO     DATA.inputs.CO];
            DATA.params_evol.MUx_EVOL  = [DATA.params_evol.MUx_EVOL  DATA.inputs.MUx    DATA.inputs.MUx];
            DATA.params_evol.MUy_EVOL  = [DATA.params_evol.MUy_EVOL  DATA.inputs.MUy    DATA.inputs.MUy];
            DATA.params_evol.MUz_EVOL  = [DATA.params_evol.MUz_EVOL  DATA.inputs.MUz    DATA.inputs.MUz];
            DATA.params_evol.K_N_EVOL  = [DATA.params_evol.K_N_EVOL  DATA.inputs.K_N    DATA.inputs.K_N];
            DATA.params_evol.K_T_EVOL  = [DATA.params_evol.K_T_EVOL  DATA.inputs.K_T    DATA.inputs.K_T];
            DATA.params_evol.L_EVOL    = [DATA.params_evol.L_EVOL    DATA.inputs.L      DATA.inputs.L];
            DATA.params_evol.DT_EVOL   = [DATA.params_evol.DT_EVOL   DATA.inputs.DT_SIM DATA.inputs.DT_SIM];

            ii = ii + 1;
            JOBFOUND = JOBFOUND + 1;
            LASTJOB  = JOBNUM;
        end
    elseif (JOBNUM > JOBNUMMAX)
        CONTINUE = 0;
        disp(['found ',num2str(JOBFOUND),' results']);
    end
    JOBNUM   = JOBNUM + 1;
end

if(JOBFOUND == 0)
    disp('no results found, please verify the paths');
    return;
else
    %% Build grids

    Nky = numel(ky); 
    if Nky > 1
        dky = ky(2); 
        Ly = 2*pi/dky;   
    else
        dky = 0;
        Ly  = 0;   
    end
    [~,iky0] = min(abs(ky)); 
    Ny = 2*Nky-1;  
    y  = linspace(-Ly/2,Ly/2,Ny+1); y = y(1:end-1);

    Nkx = numel(kx);
    if Nkx > 1
        dkx = kx(2);
        Lx = 2*pi/dkx;
    else
        dkx = 0;
        Lx  = 0;
    end    
    [~,ikx0] = min(abs(kx));
    Nx = Nkx;      
    x  = linspace(-Lx/2,Lx/2,Nx+1); x = x(1:end-1);

    Nz = numel(z);

    [KX,KY] = meshgrid(kx,ky);
    KPERP2 = KX.^2+KY.^2;
    %% Add everything in output structure
    % scaling
    DATA.scale = 1;%(1/Nx/Ny)^2;
    % Fields
    DATA.GGAMMA_RI = GGAMMA_; DATA.PGAMMA_RI = PGAMMA_; DATA.HFLUX_X = HFLUX_;
    if(KIN_E)
    DATA.HFLUX_XE = HFLUXE_;
    end
    DATA.Ts0D = Ts0D_;
    DATA.KIN_E=KIN_E;
    % grids
    DATA.grids.Parray = Parray;
    DATA.grids.Jarray = Jarray;
    DATA.grids.kx = kx; DATA.grids.ky = ky; DATA.grids.z = z; DATA.grids.Npol = -z(1)/pi;
    DATA.grids.x  = x;  DATA.grids.y  = y;
    DATA.grids.ikx0 = ikx0; DATA.grids.iky0 = iky0;
    DATA.grids.Nx = Nx; DATA.grids.Ny = Ny; DATA.grids.Nz = Nz; DATA.grids.Nkx = Nkx; DATA.grids.Nky = Nky; 
    DATA.grids.Np = numel(Parray);DATA.grids.Nj = numel(Jarray);
    if(numel(Parray)>1)
        DATA.grids.deltap = Parray(2)-Parray(1);
    else
        DATA.grids.deltap = 1;
    end
    DATA.dir      = DIRECTORY;
    DATA.localdir = DIRECTORY;
    DATA.param_title=['$\nu_{',DATA.inputs.CONAME,'}=$', num2str(DATA.inputs.NU), ...
        ', $\kappa_{Ni}=$',num2str(DATA.inputs.K_N),', $\kappa_{Ti}=$',num2str(DATA.inputs.K_T),...
        ', $L=',num2str(DATA.inputs.L),'$, $N=',...
        num2str(DATA.grids.Nx),'$, $(P,J)=(',num2str(DATA.grids.Parray(end)),',',...
        num2str(DATA.grids.Jarray(end)),')$,',' $\mu_{hd}=$(',num2str(DATA.inputs.MUx),...
        ',',num2str(DATA.inputs.MUy),')'];
    DATA.paramshort = [num2str(DATA.grids.Np),'x',num2str(DATA.grids.Nj),'x',...
        num2str(DATA.grids.Nkx),'x',num2str(DATA.grids.Nky),'x',num2str(DATA.grids.Nz)];
    JOBNUM = LASTJOB;

    filename = sprintf([DIRECTORY,'outputs_%.2d.h5'],JOBNUM);
end
end