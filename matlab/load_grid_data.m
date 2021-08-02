function [ pe, je, pi, ji, kx, ky, z ] = load_grid_data_3D( filename )
%LOAD_GRID_DATA stored in a hdf5 result file from HeLaZ
    pe    = h5read(filename,'/data/grid/coordp_e');
    je    = h5read(filename,'/data/grid/coordj_e');
    pi    = h5read(filename,'/data/grid/coordp_i');
    ji    = h5read(filename,'/data/grid/coordj_i');
    kx    = h5read(filename,'/data/grid/coordkx');
    ky    = h5read(filename,'/data/grid/coordky');
    z     = h5read(filename,'/data/grid/coordz');
end