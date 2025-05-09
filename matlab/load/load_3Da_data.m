function [ data, time, dt ] = load_3Da_data( filename, variablename )
%LOAD_3D_DATA load a 3D variable stored in a hdf5 result file from HeLaZ
    time  = h5read(filename,'/data/var3d/time');
    dt    = h5readatt(filename,'/data/input/basic','dt');
    cstart= h5readatt(filename,'/data/input/basic','start_iframe3d'); 
    Na    = h5readatt(filename,'/data/input/model','Na'); 
    Np    = h5readatt(filename,'/data/input/grid', 'Np'); 
    Nj    = h5readatt(filename,'/data/input/grid', 'Nj'); 
    Nky   = h5readatt(filename,'/data/input/grid', 'Nky'); 
    Nkx   = h5readatt(filename,'/data/input/grid', 'Nkx'); 
    Nz    = h5readatt(filename,'/data/input/grid', 'Nz'); 

    
    % Find array size by loading the first output
    tmp   = h5read(filename,['/data/var3d/',variablename,'/', num2str(cstart+1,'%06d')]);
    try % check if it is complex or real
        sz  = size(tmp.real);
        cmpx = 1;
    catch
        sz  = size(tmp);
        cmpx = 0;
    end
    % add a dimension if nz=1 or na=1
    % if Na == 1
    %     sz = [1 sz];
    % end
    if Nz == 1
        sz = [sz 1];
    end
    if Np == 1
        sz = [sz 1];
    end
    % add time dimension
    data     = zeros([sz numel(time)]);
    sz_t   = size(data(:,:,:,:,1));
    for it = 1:numel(time)
        tmp         = h5read(filename,['/data/var3d/',variablename,'/', num2str(cstart+it,'%06d')]);
        if cmpx
            data(:,:,:,:,it) = reshape(tmp.real,sz_t) + 1i * reshape(tmp.imaginary,sz_t);
        else
            data(:,:,:,:,it) = reshape(tmp,sz_t);
        end
    end

end

