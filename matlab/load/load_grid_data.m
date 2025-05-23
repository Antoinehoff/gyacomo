function [ p, j, kx, ky, z ] = load_grid_data( filename )
%LOAD_GRID_DATA stored in a hdf5 result file from HeLaZ
    p     = h5read(filename,'/data/grid/coordp');
    j     = h5read(filename,'/data/grid/coordj');
    kx    = h5read(filename,'/data/grid/coordkx');
    ky    = h5read(filename,'/data/grid/coordky');
    z     = h5read(filename,'/data/grid/coordz');
    if (numel(z) == 1) 
        z = 0;
    end
end