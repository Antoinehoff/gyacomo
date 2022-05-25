system('cd ../results/debug/test_diag; mpirun -np 8 ./helaz3_dbg 2 2 2; cd $HOME/HeLaZ/wk');


filename = '../results/debug/test_diag/outputs_00.h5';

f_        = h5read(filename,'/data/var3d/phi/000004/');
f_gatherv = h5read(filename,'/data/var3d/phi_gatherv/000004/');
f_        = f_.real        + 1i*f_.imaginary;
f_gatherv = f_gatherv.real + 1i*f_gatherv.imaginary;
err = sum(sum(sum(abs(f_-f_gatherv))))/sum(sum(sum(abs(f_))));
disp(['error_phi = ',sprintf('%2.2e',err)]);    

f_        = h5read(filename,'/data/var3d/Ni00/000004/');
f_gatherv = h5read(filename,'/data/var3d/Ni00_gatherv/000004/');
f_        = f_.real        + 1i*f_.imaginary;
f_gatherv = f_gatherv.real + 1i*f_gatherv.imaginary;
err = sum(sum(sum(abs(f_-f_gatherv))))/sum(sum(sum(abs(f_))));
disp(['error_Ni00 = ',sprintf('%2.2e',err)]);    
