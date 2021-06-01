% Write the input script "fort.90" with desired parameters
INPUT = 'setup_and_run.sh';
fid = fopen(INPUT,'wt');

fprintf(fid,[...
'#!/bin/bash\n',...
'mkdir -p $CINECA_SCRATCH/HeLaZ/wk\n',...
...
'cd $CINECA_SCRATCH/HeLaZ/wk/\n',...
...
'mkdir -p ', BASIC.RESDIR,'\n',...
'cd ',BASIC.RESDIR,'\n',...
'cp $HOME/HeLaZ/wk/fort.90 .\n',...
'cp $HOME/HeLaZ/wk/batch_script.sh .\n',...
...
'sbatch batch_script.sh\n',...
'echo tail -f $CINECA_SCRATCH/HeLaZ',BASIC.RESDIR(3:end),'out.txt']);

fclose(fid);
system(['cp setup_and_run.sh ',BASIC.RESDIR,'/.']);

% Write the sbatch script
INPUT = 'batch_script.sh';
fid = fopen(INPUT,'wt');

fprintf(fid,[...
'#!/bin/bash\n',...
'#SBATCH --job-name=',CLUSTER.JNAME,'\n',...
'#SBATCH --time=', CLUSTER.TIME,'\n',...
'#SBATCH --nodes=', CLUSTER.NODES,'\n',...
'#SBATCH --cpus-per-task=', CLUSTER.CPUPT,'\n',...
'#SBATCH --ntasks-per-node=', CLUSTER.NTPN,'\n',...
'#SBATCH --mem=', CLUSTER.MEM,'\n',...
'#SBATCH --error=err.txt\n',...
'#SBATCH --output=out.txt\n',...
'#SBATCH --account=FUA35_TSVVT421\n',...
'#SBATCH --partition=skl_fua_',CLUSTER.PART,'\n',...
'module load autoload hdf5 fftw\n',...
'srun --cpu-bind=cores ./../../../bin/',EXECNAME,' ',num2str(NP_P),' ',num2str(NP_KR)]);

fclose(fid);
system(['cp batch_script.sh ',BASIC.RESDIR,'/.']);

system('scp {fort.90,setup_and_run.sh,batch_script.sh} ahoffman@login.marconi.cineca.it:/marconi/home/userexternal/ahoffman/HeLaZ/wk > trash.txt');
system('rm trash.txt');