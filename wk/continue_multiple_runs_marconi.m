%% Paste the list of continue_run calls
continue_run('/marconi_scratch/userexternal/ahoffman/HeLaZ/results/v2.5_P_10_J_5/200x100_L_120_P_10_J_5_eta_0.5_nu_1e+00_SGGK_CLOS_0_mu_2e-02/out.txt')
continue_run('/marconi_scratch/userexternal/ahoffman/HeLaZ/results/v2.5_P_10_J_5/200x100_L_120_P_10_J_5_eta_0.6_nu_1e+00_SGGK_CLOS_0_mu_2e-02/out.txt')
continue_run('/marconi_scratch/userexternal/ahoffman/HeLaZ/results/v2.5_P_10_J_5/200x100_L_120_P_10_J_5_eta_0.7_nu_1e+00_SGGK_CLOS_0_mu_2e-02/out.txt')
continue_run('/marconi_scratch/userexternal/ahoffman/HeLaZ/results/v2.5_P_10_J_5/200x100_L_120_P_10_J_5_eta_0.8_nu_1e+00_SGGK_CLOS_0_mu_2e-02/out.txt')
continue_run('/marconi_scratch/userexternal/ahoffman/HeLaZ/results/v2.5_P_10_J_5/200x100_L_120_P_10_J_5_eta_0.6_nu_5e-01_DGGK_CLOS_0_mu_2e-02/out.txt')
continue_run('/marconi_scratch/userexternal/ahoffman/HeLaZ/results/v2.5_P_10_J_5/200x100_L_120_P_10_J_5_eta_0.7_nu_5e-01_DGGK_CLOS_0_mu_2e-02/out.txt')
continue_run('/marconi_scratch/userexternal/ahoffman/HeLaZ/results/v2.5_P_10_J_5/200x100_L_120_P_10_J_5_eta_0.8_nu_5e-01_DGGK_CLOS_0_mu_2e-02/out.txt')
continue_run('/marconi_scratch/userexternal/ahoffman/HeLaZ/results/v2.5_P_10_J_5/200x100_L_120_P_10_J_5_eta_0.9_nu_5e-01_DGGK_CLOS_0_mu_2e-02/out.txt')
continue_run('/marconi_scratch/userexternal/ahoffman/HeLaZ/results/v2.5_P_10_J_5/200x100_L_120_P_10_J_5_eta_0.9_nu_1e-01_DGGK_CLOS_0_mu_2e-02/out.txt')

%% Functions to modify preexisting fort.90 input file and launch on marconi
function [] = continue_run(outfilename)
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
        line = A{33};
        line = line(end-2:end);
        if(line(1) == '='); line = line(end); end;
        J2L = str2num(line) + 1;
    end
    % Change job 2 load in fort.90
    A{33} = ['  job2load      = ',num2str(J2L)];
    disp(A{33})
    
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