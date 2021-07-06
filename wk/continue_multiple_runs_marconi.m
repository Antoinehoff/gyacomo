%% Paste the list of continue_run calls
continue_run('/marconi_scratch/userexternal/ahoffman/HeLaZ/results/v2.7_P_6_J_3/200x100_L_60_P_6_J_3_eta_0.6_nu_1e-03_SGGK_CLOS_0_mu_1e-03/out.txt')

%% Functions to modify preexisting fort.90 input file and launch on marconi
function [] = continue_run(outfilename)
    EXECNAME = 'helaz_2.73';
    %% CLUSTER PARAMETERS
    CLUSTER.PART  = 'prod';     % dbg or prod
    CLUSTER.TIME  = '24:00:00'; % allocation time hh:mm:ss
    if(strcmp(CLUSTER.PART,'dbg')); CLUSTER.TIME  = '00:30:00'; end;
    CLUSTER.MEM   = '64GB';     % Memory
    CLUSTER.JNAME = 'HeLaZ';% Job name
    NP_P          = 2;          % MPI processes along p  
    NP_KR         = 24;         % MPI processes along kr
    % Compute processes distribution
    Ntot = NP_P * NP_KR;
    Nnodes = ceil(Ntot/48);
    Nppn   = Ntot/Nnodes; 
    CLUSTER.NODES =  num2str(Nnodes);  % MPI process along p
    CLUSTER.NTPN  =  num2str(Nppn); % MPI process along kr
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
    if( numel(A{5}) ==  numel('  RESTART   = .false.') )
        A{5} = sprintf('  RESTART   = .true.');
        J2L = 0;
    else
        line = A{35};
        line = line(end-2:end);
        if(line(1) == '='); line = line(end); end;
        J2L = str2num(line);
    end
    % Change job 2 load in fort.90
    A{35} = ['  job2load      = ',num2str(0)];
    disp(A{35})
    % Change time step
    A{3} = ['  dt     = 0.01'];
    % Increase endtime
    A{4} = ['  tmax      = 10000'];
    % Put non linear term back
    A{41} = ['  NL_CLOS = -1'];
    % change HD
    line_= A{43};
    mu_old = str2num(line_(13:end));
    A{43} = ['  mu      = ',num2str(1e-3)];

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