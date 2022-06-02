% system('cd ../results/dev/test_diag; ./helaz3_dbg; cd $HOME/HeLaZ/wk');
% system('cd ../results/dev/test_diag; mpirun -np 2 ./helaz3_dbg 1 2 1; cd $HOME/HeLaZ/wk');
system('cd ../results/dev/test_diag; mpirun -np 8 ./helaz3_dbg 2 2 2; cd $HOME/HeLaZ/wk');


filename = '../results/dev/test_diag/outputs_00.h5';

% test phi
f_        = h5read(filename,'/data/var3d/phi/000004/');
f_gatherv = h5read(filename,'/data/var3d/phi_gatherv/000004/');
f_        = f_.real        + 1i*f_.imaginary;
f_gatherv = f_gatherv.real + 1i*f_gatherv.imaginary;
err = sum(sum(sum(abs(f_-f_gatherv))))/sum(sum(sum(abs(f_))));
disp(['error_phi = ',sprintf('%2.2e',err)]);    

% test Ni00
f_        = h5read(filename,'/data/var3d/Ni00/000004/');
f_gatherv = h5read(filename,'/data/var3d/Ni00_gatherv/000004/');
f_        = f_.real        + 1i*f_.imaginary;
f_gatherv = f_gatherv.real + 1i*f_gatherv.imaginary;
err = sum(sum(sum(abs(f_-f_gatherv))))/sum(sum(sum(abs(f_))));
disp(['error_Ni00 = ',sprintf('%2.2e',err)]);    

% test Nipj
f_        = h5read(filename,'/data/var5d/moments_i/000004/');
f_gatherv = h5read(filename,'/data/var5d/moments_i_gatherv/000004/');
f_        = f_.real        + 1i*f_.imaginary;
f_gatherv = f_gatherv.real + 1i*f_gatherv.imaginary;
err = sum(sum(sum(sum(sum(abs(f_-f_gatherv))))))/sum(sum(sum(sum(sum(abs(f_))))));
disp(['error_Nipj = ',sprintf('%2.2e',err)]);    
