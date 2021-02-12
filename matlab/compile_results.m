CONTINUE = 1;
JOBNUM   = 0; JOBFOUND = 0;
Nipj_    = []; Nepj_    = [];
Ni00_    = []; Ne00_    = [];
PHI_     = [];
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
        sz = size(Nepj); Pe_new= sz(1); Je_new= sz(2);
        sz = size(Nipj); Pi_new= sz(1); Ji_new= sz(2);
        % If a degree is larger than previous job, put them in a larger array
        if (sum([Pe_new, Je_new, Pi_new, Ji_new]>[Pe_old, Je_old, Pi_old, Ji_old]) >= 1)
            tmp = Nipj_; sz = size(tmp);
            Nipj_ = zeros(cat(1,[Pi_new,Ji_new]',sz(3:end)')');
            Nipj_(1:Pi_old,1:Ji_old,:,:,:) = tmp;
            tmp = Nepj_; sz = size(tmp);
            Nepj_ = zeros(cat(1,[Pe_new,Je_new]',sz(3:end)')');
            Nepj_(1:Pe_old,1:Je_old,:,:,:) = tmp;
            tmp = Sipj_; sz = size(tmp);
            Sipj_ = zeros(cat(1,[Pi_new,Ji_new]',sz(3:end)')');
            Sipj_(1:Pi_old,1:Ji_old,:,:,:) = tmp;
            tmp = Sepj_; sz = size(tmp);
            Sepj_ = zeros(cat(1,[Pe_new,Je_new]',sz(3:end)')');
            Sepj_(1:Pe_old,1:Je_old,:,:,:) = tmp;
        end
        
        
        Nipj_ = cat(5,Nipj_,Nipj);
        Nepj_ = cat(5,Nepj_,Nepj);
        Ni00_ = cat(3,Ni00_,Ni00);
        Ne00_ = cat(3,Ne00_,Ne00);
        PHI_  = cat(3,PHI_,PHI);
        Ts2D_   = cat(1,Ts2D_,Ts2D);
        Ts5D_   = cat(1,Ts5D_,Ts5D);
        
        Sipj_ = cat(5,Sipj_,Sipj);
        Sepj_ = cat(5,Sepj_,Sepj);

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
Nipj = Nipj_; Nepj = Nepj_; Ts5D = Ts5D_;
Ni00 = Ni00_; Ne00 = Ne00_; PHI = PHI_; Ts2D = Ts2D_;
clear Nipj_ Nepj_ Ni00_ Ne00_ PHI_ Ts2D_ Ts5D_

Sipj = Sipj_; Sepj = Sepj_;
clear Sipj_ Sepj_
JOBNUM = LASTJOB
filename = sprintf([BASIC.RESDIR,'outputs_%.2d.h5'],JOBNUM);