function [ data, time, dt ] = load_5D_data( filename, variablename )
%LOAD_5D_DATA load a 5D variable stored in a hdf5 result file from HeLaZ
    time  = h5read(filename,'/data/var5d/time');
    if strcmp(variablename,'moments_e') || strcmp(variablename,'Sepj')
        p     = h5read(filename,'/data/grid/coordp_e');
        j     = h5read(filename,'/data/grid/coordj_e');
    else
        p     = h5read(filename,'/data/grid/coordp_i');
        j     = h5read(filename,'/data/grid/coordj_i');
    end
    kx    = h5read(filename,'/data/grid/coordkx');
    ky    = h5read(filename,'/data/grid/coordky');
    z     = h5read(filename,'/data/grid/coordz');

    dt    = h5readatt(filename,'/data/input','dt');
    cstart= h5readatt(filename,'/data/input','start_iframe5d'); 
    
    data  = zeros(numel(p),numel(j),numel(ky),numel(kx),numel(z),numel(time));
    
    for it = 1:numel(time)
        tmp          = h5read(filename,['/data/var5d/', variablename,'/', num2str(cstart+it,'%06d')]);
        data(:,:,:,:,:,it) = tmp.real + 1i * tmp.imaginary;
    end
end