function [ data, time, dt ] = load_2Da_data( filename, variablename )
%LOAD_2Da_DATA load a 2D variable stored in a hdf5 result file
    time  = h5read(filename,'/data/var2d/time');
    dt    = h5readatt(filename,'/data/input/basic','dt');
    cstart= h5readatt(filename,'/data/input/basic','start_iframe2d'); 
    Na    = h5readatt(filename,'/data/input/model','Na'); 
    Np    = h5readatt(filename,'/data/input/grid', 'Np'); 
    Nj    = h5readatt(filename,'/data/input/grid', 'Nj'); 
    Nky   = h5readatt(filename,'/data/input/grid', 'Nky'); 
    Nkx   = h5readatt(filename,'/data/input/grid', 'Nkx'); 
    Nz    = h5readatt(filename,'/data/input/grid', 'Nz'); 

    
    % Find array size by loading the first output
    tmp   = h5read(filename,['/data/var2d/',variablename,'/', num2str(cstart+1,'%06d')]);
    try % check if it is complex or real
        sz  = size(tmp.real);
        cmpx = 1;
    catch
        sz  = size(tmp);
        cmpx = 0;
    end
    % add a dimension if nz=1 or na=1
%     if Na == 1
%         sz = [1 sz];
%     end
    if Nz == 1
        sz = [sz 1];
    end
    % add time dimension
    data     = zeros([sz numel(time)]);
    
    for it = 1:numel(time)
        tmp         = h5read(filename,['/data/var2d/',variablename,'/', num2str(cstart+it,'%06d')]);
        if cmpx
            data(:,:,:,it) = reshape(tmp.real,sz) + 1i * reshape(tmp.imaginary,sz);
        else
            data(:,:,:,it) = reshape(tmp,sz);
        end
    end

end

