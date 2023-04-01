function [ data, time, dt ] = load_3D_data( filename, variablename )
%LOAD_3D_DATA load a 3D variable stored in a hdf5 result file from HeLaZ
    time     = h5read(filename,'/data/var3d/time');
    dt    = h5readatt(filename,'/data/input/basic','dt');
    cstart= h5readatt(filename,'/data/input/basic','start_iframe3d'); 
    
    % Find array size by loading the first output
    tmp   = h5read(filename,['/data/var3d/',variablename,'/', num2str(cstart+1,'%06d')]);
    try % check if it is complex or real
        sz  = size(tmp.real);
        cmpx = 1;
    catch
        sz  = size(tmp);
        cmpx = 0;
    end
    % add a z dimension even if 2D
    if(numel(sz) == 2)
        sz = [sz 1];
    end
    % add time dimension
    sz    = [sz numel(time)];
    data     = zeros(sz);
    
    for it = 1:numel(time)
        tmp         = h5read(filename,['/data/var3d/',variablename,'/', num2str(cstart+it,'%06d')]);
        if cmpx
            switch numel(sz)
                case(3)
                data(:,:,1,it) = tmp.real + 1i * tmp.imaginary;
                case(4)
                data(:,:,:,it) = tmp.real + 1i * tmp.imaginary;
                case(5)
                data(:,:,:,:,it) = tmp.real + 1i * tmp.imaginary;
            end
        else
            switch numel(sz)
                case(3)
                data(:,:,1,it) = tmp;
                case(4)
                data(:,:,:,it) = tmp;
                case(5)
                data(:,:,:,:,it) = tmp;
            end
        end
    end

end

