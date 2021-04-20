% Write the input script "fort.90" with desired parameters
INPUT = 'setup_and_run.sh';
fid = fopen(INPUT,'wt');

fprintf(fid,[...
'#!/bin/bash\n',...
'mkdir -p $SCRATCH/HeLaZ/wk\n',...
...
'cd $SCRATCH/HeLaZ/wk/\n',...
...
'mkdir -p ', BASIC.RESDIR,'\n',...
'cd ',BASIC.RESDIR,'\n',...
'cp $HOME/HeLaZ/wk/fort.90 .\n',...
'cp $HOME/HeLaZ/wk/batch_script.sh .\n',...
...
'jid=$(sbatch batch_script.sh)\n',...
'echo $jid\n',...
'echo to check output log :\n',...
'echo tail -f $SCRATCH/HeLaZ/results/',BASIC.SIMID,'/',BASIC.PARAMS,'/out.txt']);

fclose(fid);
system(['cp setup_and_run.sh ',BASIC.RESDIR,'/.']);

% Write the sbatch script
INPUT = 'batch_script.sh';
fid = fopen(INPUT,'wt');

fprintf(fid,[...
'#!/bin/bash\n',...
'#SBATCH --job-name="',CLUSTER.JNAME,'"\n',...
'#SBATCH --time=', CLUSTER.TIME,'\n',...
'#SBATCH --nodes=', CLUSTER.NODES,'\n',...
'#SBATCH --cpus-per-task=', CLUSTER.CPUPT,'\n',...
'#SBATCH --ntasks-per-node=', CLUSTER.NTPN,'\n',...
'#SBATCH --ntasks-per-core=', CLUSTER.NTPC,'\n',...
'#SBATCH --mem=', CLUSTER.MEM,'\n',...
'#SBATCH --error=err.txt\n',...
'#SBATCH --output=out.txt\n',...
'#SBATCH --account="s882"\n',...
'#SBATCH --constraint=mc\n',...
'#SBATCH --hint=nomultithread\n',...
'#SBATCH --partition=',CLUSTER.PART,'\n',...
'export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n',...
...% '#SBATCH --job-name=',PARAMS,'\n\n',...
'module purge\n',...
'module load PrgEnv-intel\n',...
'module load cray-hdf5-parallel\n',...
'module load cray-mpich\n',...
'module load craype-x86-skylake\n',...
'module load cray-fftw\n',...
'srun --cpu-bind=cores ./../../../bin/helaz ',num2str(NP_P),' ',num2str(NP_KR)]);
%'srun ./../../../bin/helaz']);

fclose(fid);
system(['cp batch_script.sh ',BASIC.RESDIR,'/.']);

system('scp {fort.90,setup_and_run.sh,batch_script.sh} ahoffman@ela.cscs.ch:HeLaZ/wk');