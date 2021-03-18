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
    kr    = h5read(filename,'/data/grid/coordkr');
    kz    = h5read(filename,'/data/grid/coordkz');

    dt    = h5readatt(filename,'/data/input','dt');
    cstart= h5readatt(filename,'/data/input','start_iframe5d'); 
    
    data  = zeros(numel(p),numel(j),numel(kr),numel(kz),numel(time));
    
    for it = 1:numel(time)
        tmp          = h5read(filename,['/data/var5d/', variablename,'/', num2str(cstart+it,'%06d')]);
        data(:,:,:,:,it) = tmp.real + 1i * tmp.imaginary;
    end
end