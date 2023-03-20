function [ data, time, dt ] = load_5D_data( filename, variablename )
%LOAD_5D_DATA load a 5D variable stored in a hdf5 result file from HeLaZ
    time  = h5read(filename,'/data/var5d/time');
    na    = h5readatt(filename,'/data/input/model','Na');
    p     = h5read(filename,'/data/grid/coordp');
    j     = h5read(filename,'/data/grid/coordj');
    kx    = h5read(filename,'/data/grid/coordkx');
    ky    = h5read(filename,'/data/grid/coordky');
    z     = h5read(filename,'/data/grid/coordz');

    dt    = h5readatt(filename,'/data/input/basic','dt');
    cstart= h5readatt(filename,'/data/input/basic','start_iframe5d'); 
    
    data  = zeros(na,numel(p),numel(j),numel(ky),numel(kx),numel(z),numel(time));
    
    for it = 1:numel(time)
        tmp          = h5read(filename,['/data/var5d/', variablename,'/', num2str(cstart+it,'%06d')]);
        data(:,:,:,:,:,:,it) = tmp.real + 1i * tmp.imaginary;
    end
end