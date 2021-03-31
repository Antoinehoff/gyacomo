INPUT = [COSOlver_path,'fort.90'];
fid = fopen(INPUT,'wt');

fprintf(fid,'&BASIC\n');
fprintf(fid,['pmaxe = ', num2str(COSOLVER.pmaxe),'\n']);
fprintf(fid,['jmaxe = ', num2str(COSOLVER.jmaxe),'\n']);
fprintf(fid,['pmaxi = ', num2str(COSOLVER.pmaxi),'\n']);
fprintf(fid,['jmaxi = ', num2str(COSOLVER.jmaxi),'\n']);
fprintf(fid,'JEmaxx=12\n');
fprintf(fid,'PMmaxx=16\n');
fprintf(fid,'\n');

fprintf(fid,['neFLR=',num2str(COSOLVER.neFLR),'\n']);
fprintf(fid,['neFLRs =',num2str(COSOLVER.neFLRs),'\n']);
fprintf(fid,['npeFLR=',num2str(COSOLVER.npeFLR),'\n']);
fprintf(fid,['niFLR=',num2str(COSOLVER.niFLR),'\n']);
fprintf(fid,['niFLRs=',num2str(COSOLVER.niFLRs),'\n']);
fprintf(fid,['npiFLR=',num2str(COSOLVER.npiFLR),'\n']);
fprintf(fid,'\n');
fprintf(fid,['eecolls=.true.','\n']);
fprintf(fid,['iicolls=.true.','\n']);
fprintf(fid,['eicolls=.true.','\n']);
fprintf(fid,['iecolls=.true.','\n']);
fprintf(fid,'/\n');
fprintf(fid,'\n');

fprintf(fid,'&BASIS_TRANSFORMATION_PAR\n');
fprintf(fid,['T5dir = ','''/misc/coeffs_backup/T5src/''','\n']);
fprintf(fid,['T4dir = ','''/misc/T4/NNT4_L000x200_K000x200_P000x200_J000x200/''','\n']);
fprintf(fid,'idxT4max = 30\n');
fprintf(fid,'idxT5max = 0\n');
fprintf(fid,'IFT4 = .true.\n');
fprintf(fid,'IFT5 = .false.\n');
fprintf(fid,'/\n');
fprintf(fid,'\n');

fprintf(fid,'&MODEL_PAR \n');
fprintf(fid,'nu=1\n');
fprintf(fid,'mu=0.023338\n');
fprintf(fid,'tau=1\n');
fprintf(fid,['kperp=',num2str(COSOLVER.kperp,16),'\n']);
fprintf(fid,'/\n');
fprintf(fid,'\n');

fprintf(fid,'&OPERATOR_MODEL \n');
fprintf(fid,['ETEST=', num2str(COSOLVER.co),'\n']);
fprintf(fid,['EBACK=', num2str(COSOLVER.co),'\n']);
fprintf(fid,['ITEST=', num2str(COSOLVER.co),'\n']);
fprintf(fid,['IBACK=', num2str(COSOLVER.co),'\n']);
fprintf(fid,'\n');
fprintf(fid,['ESELF=', num2str(COSOLVER.co),'\n']);
fprintf(fid,['ISELF=', num2str(COSOLVER.co),'\n']);
fprintf(fid,'\n');
fprintf(fid,['GKE = ', num2str(COSOLVER.gk),'\n']);
fprintf(fid,['GKI = ', num2str(COSOLVER.gk),'\n']);
fprintf(fid,'DKTEST = .F.\n');
fprintf(fid,'DKBACK = .F.\n');
fprintf(fid,'ADDTEST = .T.\n');
fprintf(fid,'ADDBACK = .T.\n');
fprintf(fid,'only_gk_part = .false.\n');
fprintf(fid,'only_symmetric = .false.\n');
fprintf(fid,'/\n');
fprintf(fid,'\n');

fprintf(fid,'&MPI\n');
fprintf(fid,'Nprocs_j_ii = 1\n');
fprintf(fid,'MPI_balanced = .false.\n');
fprintf(fid,'/\n');
fprintf(fid,'\n');

fprintf(fid,'&OUTPUT_PAR\n');
fprintf(fid,'OVERWRITE=.true.\n');
fprintf(fid,'nsave_ei=0\n');
fprintf(fid,'nsave_ie=0\n');
fprintf(fid,'nsave_ee=0\n');
fprintf(fid,'nsave_ii=2\n');
fprintf(fid,'ifrestart=.false.\n');
fprintf(fid,['suffix_resfile=','''''','\n']);
fprintf(fid,['outdir = ','''../../../HeLaZ/iCa''','\n']);
fprintf(fid,'/\n');
fprintf(fid,'\n');


fclose(fid);