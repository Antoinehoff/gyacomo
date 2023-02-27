function [ data, time, dt ] = load_3D_data( filename, variablename )
%LOAD_3D_DATA load a 3D variable stored in a hdf5 result file from HeLaZ
    time     = h5read(filename,'/data/var3d/time');
    dt    = h5readatt(filename,'/data/input','dt');
    cstart= h5readatt(filename,'/data/input','start_iframe3d'); 
    
    % Find array size by loading the first output
    tmp   = h5read(filename,['/data/var3d/',variablename,'/', num2str(cstart+1,'%06d')]);
    try % check if it is complex or real
        sz  = size(tmp.real);
        cmpx = 1;
    catch
        sz  = size(tmp);
        cmpx = 0;
    end
    % add time dimension
    sz    = [sz numel(time)];
    data     = zeros(sz);
    
    for it = 1:numel(time)
        tmp         = h5read(filename,['/data/var3d/',variablename,'/', num2str(cstart+it,'%06d')]);
        if cmpx
            if(numel(sz) == 3)
                data(:,:,it) = tmp.real + 1i * tmp.imaginary;
            else
                data(:,:,:,it) = tmp.real + 1i * tmp.imaginary;
            end
        else
            if(numel(sz) == 3)
                data(:,:,it) = tmp;
            else
                data(:,:,:,it) = tmp;
            end
        end
    end

end

