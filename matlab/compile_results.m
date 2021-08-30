CONTINUE = 1;
JOBNUM   = JOBNUMMIN; JOBFOUND = 0;
TJOB_SE  = []; % Start and end times of jobs
NU_EVOL  = []; % evolution of parameter nu between jobs
MU_EVOL  = []; % evolution of parameter mu between jobs
ETAN_EVOL= []; %
L_EVOL   = []; % 
DT_EVOL  = []; %
% FIELDS
Nipj_    = []; Nepj_    = [];
Ni00_    = []; Ne00_    = [];
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
    filename = sprintf([BASIC.MISCDIR,'outputs_%.2d.h5'],JOBNUM);
    if (exist(filename, 'file') == 2 && JOBNUM <= JOBNUMMAX)
        % Load results of simulation #JOBNUM
        load_results
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
            GGAMMA_ = cat(1,GGAMMA_,GGAMMA_RI);
            PGAMMA_ = cat(1,PGAMMA_,PGAMMA_RI);
            Ts0D_   = cat(1,Ts0D_,Ts0D);
        end
        
        if W_PHI || W_NA00
        	Ts3D_   = cat(1,Ts3D_,Ts3D);
        end
        if W_PHI
            PHI_  = cat(4,PHI_,PHI);
        end
        if W_NA00
            Ni00_ = cat(4,Ni00_,Ni00);
            Ne00_ = cat(4,Ne00_,Ne00);
        end
        if W_DENS
            DENS_E_ = cat(4,DENS_E_,DENS_E);
            DENS_I_ = cat(4,DENS_I_,DENS_I);
        end
        if W_TEMP
            TEMP_E_ = cat(4,TEMP_E_,TEMP_E);
            TEMP_I_ = cat(4,TEMP_I_,TEMP_I);
        end
        if W_NAPJ || W_SAPJ
            Ts5D_   = cat(1,Ts5D_,Ts5D);
        end
        if W_NAPJ
            Nipj_ = cat(6,Nipj_,Nipj);
            Nepj_ = cat(6,Nepj_,Nepj);
        end
        if W_SAPJ
     	  Sipj_ = cat(6,Sipj_,Sipj);
          Sepj_ = cat(6,Sepj_,Sepj);
        end

        % Evolution of simulation parameters
        load_params
        TJOB_SE   = [TJOB_SE Ts0D(1) Ts0D(end)]; 
        NU_EVOL   = [NU_EVOL NU NU];
        MU_EVOL   = [MU_EVOL MU MU];
        ETAN_EVOL = [ETAN_EVOL ETAN ETAN];
        L_EVOL    = [L_EVOL L L];
        DT_EVOL   = [DT_EVOL DT_SIM DT_SIM];
    
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
GGAMMA_RI = GGAMMA_; PGAMMA_RI = PGAMMA_; Ts0D = Ts0D_;
Nipj = Nipj_; Nepj = Nepj_; Ts5D = Ts5D_;
Ni00 = Ni00_; Ne00 = Ne00_; PHI = PHI_; Ts3D = Ts3D_;
DENS_E = DENS_E_; DENS_I = DENS_I_; TEMP_E = TEMP_E_; TEMP_I = TEMP_I_;
clear Nipj_ Nepj_ Ni00_ Ne00_ PHI_ Ts2D_ Ts5D_ GGAMMA_ PGAMMA_ Ts0D_

Sipj = Sipj_; Sepj = Sepj_;
clear Sipj_ Sepj_
JOBNUM = LASTJOB;
filename = sprintf([BASIC.RESDIR,'outputs_%.2d.h5'],JOBNUM);