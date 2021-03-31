CONTINUE = 1;
JOBNUM   = 0; JOBFOUND = 0;
Nipj_    = []; Nepj_    = [];
Ni00_    = []; Ne00_    = [];
GGAMMA_  = [];
PGAMMA_  = [];
PHI_     = [];
Ts0D_    = [];
Ts2D_    = [];
Ts5D_    = [];
Sipj_    = []; Sepj_    = [];
Pe_old   = 1e9; Pi_old = Pe_old; Je_old = Pe_old; Ji_old = Pe_old;

while(CONTINUE) 
    filename = sprintf([BASIC.RESDIR,'outputs_%.2d.h5'],JOBNUM);
    if exist(filename, 'file') == 2
        % Load results of simulation #JOBNUM
        load_results
        % Check polynomials degrees
        Pe_new= numel(Pe); Je_new= numel(Je);
        Pi_new= numel(Pi); Ji_new= numel(Ji);
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
        end
        
        if W_GAMMA
            GGAMMA_ = cat(1,GGAMMA_,GGAMMA_RI);
            PGAMMA_ = cat(1,PGAMMA_,PGAMMA_RI);
            Ts0D_   = cat(1,Ts0D_,Ts0D);
        end
        
        if W_PHI || W_NA00
        	Ts2D_   = cat(1,Ts2D_,Ts2D);
        end
        if W_PHI
            PHI_  = cat(3,PHI_,PHI);
        end
        if W_NA00
            Ni00_ = cat(3,Ni00_,Ni00);
            Ne00_ = cat(3,Ne00_,Ne00);
        end

        if W_NAPJ || W_SAPJ
            Ts5D_   = cat(1,Ts5D_,Ts5D);
        end
        if W_NAPJ
            Nipj_ = cat(5,Nipj_,Nipj);
            Nepj_ = cat(5,Nepj_,Nepj);
        end
        if W_SAPJ
     	  Sipj_ = cat(5,Sipj_,Sipj);
          Sepj_ = cat(5,Sepj_,Sepj);
        end

        JOBFOUND = JOBFOUND + 1;
        LASTJOB  = JOBNUM;
    elseif (JOBNUM > 20)
        CONTINUE = 0;
        disp(['found ',num2str(JOBFOUND),' results']);
    end
    JOBNUM   = JOBNUM + 1;
    Pe_old = Pe_new; Je_old = Je_new;
    Pi_old = Pi_new; Ji_old = Ji_new;
end
GGAMMA_RI = GGAMMA_; PGAMMA_RI = PGAMMA_; Ts0D = Ts0D_;
Nipj = Nipj_; Nepj = Nepj_; Ts5D = Ts5D_;
Ni00 = Ni00_; Ne00 = Ne00_; PHI = PHI_; Ts2D = Ts2D_;
clear Nipj_ Nepj_ Ni00_ Ne00_ PHI_ Ts2D_ Ts5D_ GGAMMA_ PGAMMA_ Ts0D_

Sipj = Sipj_; Sepj = Sepj_;
clear Sipj_ Sepj_
JOBNUM = LASTJOB
filename = sprintf([BASIC.RESDIR,'outputs_%.2d.h5'],JOBNUM);