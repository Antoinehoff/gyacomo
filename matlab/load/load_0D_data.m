function [ data, time, dt ] = load_0D_data( filename, variablename )
%LOAD_0D_DATA load a 0D variable stored in a hdf5 result file from HeLaZ
    time  = h5read(filename,'/data/var0d/time');
    dt    = h5readatt(filename,'/data/input/basic','dt');
    Na    = h5readatt(filename,'/data/input/model','Na');
    cstart= h5readatt(filename,'/data/input/basic','start_iframe0d'); 
    data  = h5read(filename,['/data/var0d/',variablename]);
%     data  = reshape(data,[Na, numel(time)]);
end

