function [ data, time, dt ] = load_3D_data( filename, variablename )
%LOAD_3D_DATA load a 3D variable stored in a hdf5 result file from HeLaZ
    time     = h5read(filename,'/data/var3d/time');
    dt    = h5readatt(filename,'/data/input','dt');
    cstart= h5readatt(filename,'/data/input','start_iframe3d'); 
    
    % Find array size by loading the first output
    tmp   = h5read(filename,['/data/var3d/',variablename,'/', num2str(cstart+1,'%06d')]);
    sz    = size(tmp.real);
    % add time dimension
    sz    = [sz numel(time)];
    data     = zeros(sz);
    
    for it = 1:numel(time)
        tmp         = h5read(filename,['/data/var3d/',variablename,'/', num2str(cstart+it,'%06d')]);
        data(:,:,:,it) = tmp.real + 1i * tmp.imaginary;
    end

end

