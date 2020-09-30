function [ data, p, j, kr, kz, time ] = load_5D_data( filename, variablename )
%LOAD_5D_DATA load a 5D variable stored in a hdf5 result file from HeLaZ
    time  = h5read(filename,'/data/var5d/time');
    p     = h5read(filename,['/data/var5d/',variablename,'/coordp']);
    j     = h5read(filename,['/data/var5d/',variablename,'/coordj']);
    kr    = h5read(filename,['/data/var5d/',variablename,'/coordkr']);
    kz    = h5read(filename,['/data/var5d/',variablename,'/coordkz']);
    data  = zeros(numel(p),numel(j),numel(kr),numel(kz),numel(time));
    for it = 1:numel(time)
        tmp          = h5read(filename,['/data/var5d/', variablename,'/', num2str(it,'%06d')]);
        data(:,:,:,:,it) = tmp.real + 1i * tmp.imaginary;
    end
end