%% Paste the list of continue_run calls
% continue_run('/marconi_scratch/userexternal/ahoffman/HeLaZ/results/simulation_A/300x150_L_120_P_8_J_4_eta_0.6_nu_1e-01_SGGK_mu_0e+00/out.txt')
continue_run('/marconi_scratch/userexternal/ahoffman/HeLaZ/results/simulation_B/300x150_L_120_P_8_J_4_eta_0.6_nu_5e-01_SGGK_mu_0e+00/out.txt')

%% Functions to modify preexisting fort.90 input file and launch on marconi
function [] = continue_run(outfilename)
    EXECNAME = 'helaz_3.81';
    SBATCH_CMD = 'sbatch batch_script.sh\n';
    %% CLUSTER PARAMETERS
    CLUSTER.PART  = 'prod';     % dbg or prod
    CLUSTER.TIME  = '24:00:00'; % allocation time hh:mm:ss
    if(strcmp(CLUSTER.PART,'dbg')); CLUSTER.TIME  = '00:30:00'; end;
    CLUSTER.MEM   = '64GB';     % Memory
    CLUSTER.JNAME = 'HeLaZ';% Job name
    NP_P          = 2;          % MPI processes along p  
    NP_KX         = 24;         % MPI processes along kx
    % Compute processes distribution
    Ntot = NP_P * NP_KX;
    Nnodes = ceil(Ntot/48);
    Nppn   = Ntot/Nnodes; 
    CLUSTER.NODES =  num2str(Nnodes);  % MPI process along p
    CLUSTER.NTPN  =  num2str(Nppn); % MPI process along kx
    CLUSTER.CPUPT = '1';        % CPU per task
    %%
    RESDIR = ['../',outfilename(46:end-8),'/'];
    BASIC.RESDIR = RESDIR;
    FORT90 = [RESDIR,'fort.90'];
  
    % Read txt into cell A
    fid = fopen(FORT90,'r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
    end
    fclose(fid);

    % Find previous job2load
    line    = A{38};
    J2L_old = str2num(line(17:end););
    A{38} = ['  job2load      = ',num2str(J2L_old+1)];
    disp(A{38})
    % Change time step
    line_= A{3};
    dt_old = str2num(line_(13:end));
    A{3} = ['  dt     = ',num2str(dt_old)];
    % Increase endtime
    A{4} = ['  tmax      = 20000'];
    % Change collision operator
    line_= A{42};
    CO_old = str2num(line_(13:end));
    A{42} = ['  CO      = ',num2str(2)];
    % Put non linear term back
    A{44} = ['  NL_CLOS = -1'];
    % change HD
    line_= A{46};
    mu_old = str2num(line_(13:end));
    A{46} = ['  mu      = ',num2str(0)];
    % change L
    line_= A{13};
    L_old = str2num(line_(12:end));
    A{13} = ['  Lx     = ',num2str(L_old)];
    A{15} = ['  Ly     = ',num2str(L_old)];
    % change eta N
    line_= A{56};
    etan_old = str2num(line_(13:end));
    A{56} = ['  eta_n   = ',num2str(etan_old)];
    % change eta B
    line_= A{58};
    etab_old = str2num(line_(13:end));
    A{58} = ['  eta_B   = ',num2str(etab_old)];
    % Rewrite fort.90
    fid = fopen('fort.90', 'w');
    for i = 1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break
        else
            fprintf(fid,'%s\n', A{i});
        end
    end
    % Copy fort.90 into marconi
    write_sbash_marconi
    % Launch the job
    system('ssh ahoffman@login.marconi.cineca.it sh HeLaZ/wk/setup_and_run.sh');
end
