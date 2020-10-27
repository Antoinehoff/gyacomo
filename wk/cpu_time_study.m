clear all;
CPUFREQ   = 3.2e9;%[Hz] cpu clock of the computer
CTIME_EXP = 425;  %[s] real time for N=64, P=2, J=1, DT = 1e-2, TMAX = 150
CCOST_EXP = 64*64/2*4*(2*1 + 2*1)*150/1e-2;

FACTOR    = CTIME_EXP/(CCOST_EXP/CPUFREQ);

%% Default parameters
TMAX    = 20;  % Maximal time unit
DT      = 5e-2;   % Time step
N       = 20;     % Frequency gridpoints (Nkr = N/2)
PMAXE   = 00;     % Highest electron Hermite polynomial degree
JMAXE   = 00;     % Highest ''       Laguerre ''
PMAXI   = 00;     % Highest ion      Hermite polynomial degree
JMAXI   = 00;     % Highest ''       Laguerre ''

%% CPUTIME VS DT
DTA = logspace(-3,-1,10);
NN_DT = []; CT_EST_DT = []; CT_REAL_DT = [];
for DT_ = DTA
    DT = DT_;   
    parameters_cpu_time;
    run;
    load_results;
    CCOST      = 4*N*N/2*((PMAXE+1)*(JMAXE+1) + (PMAXI+1)*(JMAXI+1))*TMAX/DT;
    NN_DT      = [NN_DT,TMAX/DT_]; 
    CT_EST_DT  = [CT_EST_DT,CCOST/CPUFREQ]; 
    CT_REAL_DT = [CT_REAL_DT,CPUTIME];
end
system('rm test_cputime*');
if 0
%%
figure
    plot(CT_EST_DT,CT_REAL_DT,'-x'); hold on;
    xlabel('$4\times N^2/2 \times(P_{e+1} J_{e+1}+ P_{i+1} J_{i+1})\times N_t/\Omega_{CPU}$'); 
    ylabel('$\tau_\textrm{CPU}$'); grid on;
end  
%% CPUTIME VS N
DT      = 5e-2;   % Time step
NA = [16 32 64 128 256 512];
NN_N = []; CT_EST_N = []; CT_REAL_N = [];
for N_ = NA
    N = N_;   
    parameters_cpu_time;
    run;
    load_results;
    CCOST     = FACTOR * N*N/2*4*((PMAXE+1)*(JMAXE+1) + (PMAXI+1)*(JMAXI+1))*TMAX/DT;
    NN_N      = [NN_N,TMAX/DT_]; 
    CT_EST_N  = [CT_EST_N,CCOST/CPUFREQ]; 
    CT_REAL_N = [CT_REAL_N,CPUTIME];
end
system('rm test_cputime*');
if 0
%%
figure
    plot(CT_EST_N,CT_REAL_N,'-x'); hold on;
    xlabel('$4\times N^2/2 \times(P_{e+1} J_{e+1}+ P_{i+1} J_{i+1})\times N_t/\Omega_{CPU}$'); 
    ylabel('$\tau_\textrm{CPU}$'); grid on;
end
%% CPUTIME VS P
DT   = 5e-2;   % Time step
N    = 20;
PA   = 0:20;
NN_P = []; CT_EST_P = []; CT_REAL_P = [];
for P_ = PA
    PMAXE = P_; PMAXI = P_;   
    parameters_cpu_time;
    run;
    load_results;
    CCOST     = FACTOR * N*N/2*4*((PMAXE+1)*(JMAXE+1) + (PMAXI+1)*(JMAXI+1))*TMAX/DT;
    NN_P      = [NN_P,TMAX/DT_]; 
    CT_EST_P  = [CT_EST_P,CCOST/CPUFREQ]; 
    CT_REAL_P = [CT_REAL_P,CPUTIME];
end
system('rm test_cputime*');
if 0
%%
figure
    plot(CT_EST_P,CT_REAL_P,'-x'); hold on;
    xlabel('$4\times N^2/2 \times\sum_a P_{a+1} J_{a+1}\times N_t/\Omega_{CPU}$[s]'); 
    ylabel('$\tau_\textrm{CPU}$[s]'); grid on;
end
%% CPUTIME VS J
DT   = 5e-2;   % Time step
N    = 20;
PMAXI= 0; PMAXE = 0;
JA   = 0:10;
NN_J = []; CT_EST_J = []; CT_REAL_J = [];
for J_ = JA
    JMAXE = J_; JMAXI = J_;   
    parameters_cpu_time;
    run;
    load_results;
    CCOST     = FACTOR * N*N/2*4*((PMAXE+1)*(JMAXE+1) + (PMAXI+1)*(JMAXI+1))*TMAX/DT;
    NN_J      = [NN_J,TMAX/DT_]; 
    CT_EST_J  = [CT_EST_J,CCOST/CPUFREQ]; 
    CT_REAL_J = [CT_REAL_J,CPUTIME];
end
system('rm test_cputime*');
if 0
%%
figure
    plot(CT_EST_P,CT_REAL_P,'-x'); hold on;
    xlabel('$4\times N^2/2 \times\sum_a P_{a+1} J_{a+1}\times N_t/\Omega_{CPU}$[s]'); 
    ylabel('$\tau_\textrm{CPU}$[s]'); grid on;
end
%% Overall comparison
figure
    loglog(CT_EST_DT,CT_REAL_DT,'-x','DisplayName','$\Delta t$'); hold on;
    plot(CT_EST_N, CT_REAL_N, '-x','DisplayName','$N$'); hold on;
    plot(CT_EST_P, CT_REAL_P, '-x','DisplayName','$P$'); hold on;
    plot(CT_EST_J, CT_REAL_J, '-x','DisplayName','$J$'); hold on;
    xlabel('$N^{2}(\sum_a P_{a+1} J_{a+1})^{1}N_s^{1} \Omega_{cpu}^{-1}$[s]'); 
    ylabel('$\tau_\textrm{CPU}$'); grid on; legend('show');
    
%%
% CTIME     = FACTOR * CCOST/CPUFREQ;
% HOURS     = floor(CTIME/3600);
% MINUTES   = floor((CTIME-3600*HOURS)/60);
% SECONDS   = floor(CTIME-3600*HOURS-60*MINUTES);
% 
% disp(['Computational cost ~',sprintf('%.1e',CCOST),' op.']);
% 
% TIMESTRING = 'Computational time ~';
% if HOURS > 0
%     TIMESTRING = [TIMESTRING,num2str(HOURS),'h '];
% end
% if MINUTES > 0
%     TIMESTRING = [TIMESTRING,num2str(MINUTES),'min '];
% end
% if (HOURS + MINUTES == 0) && (SECONDS > 0)
%     TIMESTRING = [TIMESTRING,num2str(SECONDS),'s '];
% end
% disp(TIMESTRING);
