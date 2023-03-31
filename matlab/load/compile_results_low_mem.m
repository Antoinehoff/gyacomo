function [DATA] = compile_results_low_mem(DATA,DIRECTORY,JOBNUMMIN,JOBNUMMAX)
CONTINUE = 1;
JOBNUM   = JOBNUMMIN; JOBFOUND = 0;
DATA.TJOB_SE  = []; % Start and end times of jobs
DATA.NU_EVOL  = []; % evolution of parameter nu between jobs
DATA.CO_EVOL  = []; % evolution of CO
DATA.MUx_EVOL  = []; % evolution of parameter mu between jobs
DATA.MUy_EVOL  = []; % evolution of parameter mu between jobs
DATA.MUz_EVOL  = []; % evolution of parameter mu between jobs
DATA.K_N_EVOL = []; %
DATA.K_T_EVOL = []; %
DATA.L_EVOL   = []; % 
DATA.DT_EVOL  = []; %
% Low memoery cost data
HFLUXI_   = [];
HFLUXE_   = [];
GGAMMAI_ = [];
PGAMMAI_ = [];
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
            disp(sprintf('Loading ID %.2d (%s)',JOBNUM,filename));
            % Loading from output file
            CPUTIME   = h5readatt(filename,'/data/input','cpu_time');
            DT_SIM    = h5readatt(filename,'/data/input/basic','dt');
            [Pe, Je, Pi, Ji, kx, ky, z] = load_grid_data(filename);
            W_GAMMA   = strcmp(h5readatt(filename,'/data/input/diag_par','write_gamma'),'y');
            W_HF      = strcmp(h5readatt(filename,'/data/input/diag_par','write_hf'   ),'y');
            KIN_E     = strcmp(h5readatt(filename,'/data/input/model',     'ADIAB_E' ),'n');

            if W_GAMMA
                [ GGAMMA_RI, Ts0D, ~] = load_0D_data(filename, 'gflux_xi');
                PGAMMA_RI            = load_0D_data(filename, 'pflux_xi');
                GGAMMAI_ = cat(1,GGAMMAI_,GGAMMA_RI); clear GGAMMA_RI
                PGAMMAI_ = cat(1,PGAMMAI_,PGAMMA_RI); clear PGAMMA_RI
            end

            if W_HF
                [ HFLUX_XI, Ts0D, ~] = load_0D_data(filename, 'hflux_xi');
                HFLUXI_ = cat(1,HFLUXI_,HFLUX_XI); clear HFLUX_XI
            end

            Ts0D_   = cat(1,Ts0D_,Ts0D);

            % Evolution of simulation parameters
            load_params
            DATA.TJOB_SE   = [DATA.TJOB_SE   Ts0D(1)     Ts0D(end)];
            DATA.NU_EVOL   = [DATA.NU_EVOL   DATA.NU     DATA.NU];
            DATA.CO_EVOL   = [DATA.CO_EVOL   DATA.CO     DATA.CO];
            DATA.MUx_EVOL  = [DATA.MUx_EVOL  DATA.MUx    DATA.MUx];
            DATA.MUy_EVOL  = [DATA.MUy_EVOL  DATA.MUy    DATA.MUy];
            DATA.MUz_EVOL  = [DATA.MUz_EVOL  DATA.MUz    DATA.MUz];
            DATA.K_N_EVOL  = [DATA.K_N_EVOL DATA.K_N   DATA.K_N];
            DATA.K_T_EVOL  = [DATA.K_T_EVOL DATA.K_T   DATA.K_T];
            DATA.L_EVOL    = [DATA.L_EVOL    DATA.L      DATA.L];
            DATA.DT_EVOL   = [DATA.DT_EVOL   DATA.DT_SIM DATA.DT_SIM];

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
    DATA.GGAMMA_RI = GGAMMAI_; DATA.PGAMMA_RI = PGAMMAI_; DATA.HFLUX_X = HFLUXI_;
    if(KIN_E)
    DATA.HFLUX_XE = HFLUXE_;
    end
    DATA.Ts0D = Ts0D_;
    DATA.KIN_E=KIN_E;
    % grids
    DATA.Pe = Pe; DATA.Pi = Pi; 
    DATA.Je = Je; DATA.Ji = Ji; 
    DATA.kx = kx; DATA.ky = ky; DATA.z = z; DATA.Npol = -z(1)/pi;
    DATA.x  = x;  DATA.y  = y;
    DATA.ikx0 = ikx0; DATA.iky0 = iky0;
    DATA.Nx = Nx; DATA.Ny = Ny; DATA.Nz = Nz; DATA.Nkx = Nkx; DATA.Nky = Nky; 
    DATA.Pmaxe = numel(Pe); DATA.Pmaxi = numel(Pi); DATA.Jmaxe = numel(Je); DATA.Jmaxi = numel(Ji);
    DATA.dir      = DIRECTORY;
    DATA.localdir = DIRECTORY;
    DATA.param_title=['$\nu_{',DATA.CONAME,'}=$', num2str(DATA.NU), ...
        ', $\kappa_{Ni}=$',num2str(DATA.K_N),', $\kappa_{Ti}=$',num2str(DATA.K_T),...
        ', $L=',num2str(DATA.L),'$, $N=',...
        num2str(DATA.Nx),'$, $(P,J)=(',num2str(DATA.Pi(end)),',',...
        num2str(DATA.Ji(end)),')$,',' $\mu_{hd}=$(',num2str(DATA.MUx),...
        ',',num2str(DATA.MUy),')'];
    DATA.paramshort = [num2str(DATA.Pmaxi),'x',num2str(DATA.Jmaxi),'x',...
        num2str(DATA.Nkx),'x',num2str(DATA.Nky),'x',num2str(DATA.Nz)];
    JOBNUM = LASTJOB;

    filename = sprintf([DIRECTORY,'outputs_%.2d.h5'],JOBNUM);
end
end