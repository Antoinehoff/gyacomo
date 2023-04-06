function [ data, time, dt ] = load_5Da_data( filename, variablename )
%LOAD_3D_DATA load a 3D variable stored in a hdf5 result file from HeLaZ
    time  = h5read(filename,'/data/var5d/time');
    dt    = h5readatt(filename,'/data/input/basic','dt');
    cstart= h5readatt(filename,'/data/input/basic','start_iframe3d'); 
    Na    = h5readatt(filename,'/data/input/model','Na'); 
    Np    = h5readatt(filename,'/data/input/grid', 'Np'); 
    Nj    = h5readatt(filename,'/data/input/grid', 'Nj'); 
    Nky   = h5readatt(filename,'/data/input/grid', 'Nky'); 
    Nkx   = h5readatt(filename,'/data/input/grid', 'Nkx'); 
    Nz    = h5readatt(filename,'/data/input/grid', 'Nz'); 

    
    % Find array size by loading the first output
    tmp   = h5read(filename,['/data/var5d/',variablename,'/', num2str(cstart+1,'%06d')]);
    try % check if it is complex or real
        sz  = size(tmp.real);
        cmpx = 1;
    catch
        sz  = size(tmp);
        cmpx = 0;
    end
    % add time dimension
    data     = zeros([Na,Np,Nj,Nky,Nkx,Nz,numel(time)]);
    
    for it = 1:numel(time)
        tmp         = h5read(filename,['/data/var5d/',variablename,'/', num2str(cstart+it,'%06d')]);
        if cmpx
            data(:,:,:,:,:,:,it) = reshape(tmp.real,sz) + 1i * reshape(tmp.imaginary,sz);
        else
            data(:,:,:,:,:,:,it) = reshape(tmp,sz);
        end
    end

end

