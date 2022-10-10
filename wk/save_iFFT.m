%% INPUT
simdir    = '/misc/gyacomo_outputs/results/Zpinch_rerun';
resolu    = 'UHD_512x256x2x1';
output    = 'outputs_01.h5';
filename  = [simdir,'/',resolu,'/',output]; 
fieldname = 'Ni00';

%% OUTPUT
outdir    = '.';
outname   = [fieldname,'_real_',resolu,'.h5'];
outfile   = [outdir,'/',outname];

%% Load the complex field F_c
[ F_c, Ts3D, ~] = load_3D_data(filename, fieldname);

%% Measure size and allocate the real field ff_
sz_  =size(real(fftshift(ifourier_GENE(F_c(:,:,1,1))))); % size one frame for the real dimensions
f_r  = zeros(sz_(1),sz_(2),numel(Ts3D));
sz_  = size(f_r);

%% Fill the real field f_r
for it = 1:numel(Ts3D)
    f_r(:,:,it) = squeeze(real(fftshift(ifourier_GENE(F_c(:,:,1,it)))));
end
clear F_c;

%% Save the real field
h5create(outfile,'/data',sz_,'Datatype','single')
h5write(outfile,'/data',f_r)
h5disp(outfile)
clear f_r

%% small check
test_ = h5read(outfile,'/data');
figure; imagesc(test_(:,:,end));