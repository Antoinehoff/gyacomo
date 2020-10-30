CONTINUE = 1;
JOBNUM   = 0;
Nipj_    = []; Nepj_    = [];
Ni00_    = []; Ne00_    = [];
PHI_     = [];
Ts2D_    = [];
Ts5D_    = [];

if strcmp(OUTPUTS.write_non_lin,'.true.')
Sipj_    = []; Sepj_    = [];
end

while(CONTINUE) 
    filename = sprintf([BASIC.RESDIR,'outputs_%.2d.h5'],JOBNUM);
    if exist(filename, 'file') == 2
        load_results
        Nipj_ = cat(5,Nipj_,Nipj);
        Nepj_ = cat(5,Nepj_,Nepj);
        Ni00_ = cat(3,Ni00_,Ni00);
        Ne00_ = cat(3,Ne00_,Ne00);
        PHI_  = cat(3,PHI_,PHI);
        Ts2D_   = cat(1,Ts2D_,Ts2D);
        Ts5D_   = cat(1,Ts5D_,Ts5D);
        
if strcmp(OUTPUTS.write_non_lin,'.true.')
        Sipj_ = cat(5,Sipj_,Sipj);
        Sepj_ = cat(5,Sepj_,Sepj);
end
        JOBNUM = JOBNUM + 1;
    else
        CONTINUE = 0;
        disp(['found ',num2str(JOBNUM),' results']);
    end
end
Nipj = Nipj_; Nepj = Nepj_; Ts5D = Ts5D_;
Ni00 = Ni00_; Ne00 = Ne00_; PHI = PHI_; Ts2D = Ts2D_;
clear Nipj_ Nepj_ Ni00_ Ne00_ PHI_ Ts2D_ Ts5D_

if strcmp(OUTPUTS.write_non_lin,'.true.')
Sipj = Sipj_; Sepj = Sepj_;
clear Sipj_ Sepj_
end