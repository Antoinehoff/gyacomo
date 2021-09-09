%% Paste the list of continue_run calls

continue_run('/marconi_scratch/userexternal/ahoffman/HeLaZ/results/HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_7e-02_SGGK_mu_0e+00/out.txt')

%% Functions to modify preexisting fort.90 input file and launch on marconi
function [] = continue_run(outfilename)
    EXECNAME = 'helaz_3.2';
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
    if( numel(A{5}) ==  numel('  RESTART    = .false.') )
        A{5} = sprintf('  RESTART   = .true.');
        J2L = 0;
    else
        line = A{39};
        line = line(end-2:end);
        if(line(1) == '='); line = line(end); end;
        J2L = str2num(line) + 1;
    end
    % Change job 2 load in fort.90
    A{39} = ['  job2load      = ',num2str(J2L)];
    disp(A{39})
    % Change time step
    A{3} = ['  dt     = 0.005'];
    % Increase endtime
    A{4} = ['  tmax      = 20000'];
    % Change collision operator
    line_= A{43};
    CO_old = str2num(line_(13:end));
    A{43} = ['  CO      = ',num2str(1)];
    % Put non linear term back
    A{45} = ['  NL_CLOS = -1'];
    % change HD
    line_= A{47};
    mu_old = str2num(line_(13:end));
    A{47} = ['  mu      = ',num2str(mu_old)];
    % change L
    line_= A{14};
    L_old = str2num(line_(12:end));
    A{14} = ['  Lx     = ',num2str(L_old)];
    A{16} = ['  Ly     = ',num2str(L_old)];
    % change eta N
    line_= A{57};
    etan_old = str2num(line_(13:end));
    A{57} = ['  eta_n   = ',num2str(etan_old)];
    % change eta B
    line_= A{59};
    etab_old = str2num(line_(13:end));
    A{59} = ['  eta_B   = ',num2str(etab_old)];
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
