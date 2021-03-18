function [ data, time, dt ] = load_2D_data( filename, variablename )
%LOAD_2D_DATA load a 2D variable stored in a hdf5 result file from HeLaZ
    time     = h5read(filename,'/data/var2d/time');
    kr       = h5read(filename,'/data/grid/coordkr');
    kz       = h5read(filename,'/data/grid/coordkz');
    dt    = h5readatt(filename,'/data/input','dt');
    cstart= h5readatt(filename,'/data/input','start_iframe2d'); 
    
    data     = zeros(numel(kr),numel(kz),numel(time));
    for it = 1:numel(time)
        tmp         = h5read(filename,['/data/var2d/',variablename,'/', num2str(cstart+it,'%06d')]);
        data(:,:,it) = tmp.real + 1i * tmp.imaginary;
    end

end

