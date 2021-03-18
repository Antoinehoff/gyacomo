function [ pe, je, pi, ji, kr, kz ] = load_grid_data( filename )
%LOAD_GRID_DATA stored in a hdf5 result file from HeLaZ
    pe    = h5read(filename,'/data/grid/coordp_e');
    je    = h5read(filename,'/data/grid/coordj_e');
    pi    = h5read(filename,'/data/grid/coordp_i');
    ji    = h5read(filename,'/data/grid/coordj_i');
    kr    = h5read(filename,'/data/grid/coordkr');
    kz    = h5read(filename,'/data/grid/coordkz');
end