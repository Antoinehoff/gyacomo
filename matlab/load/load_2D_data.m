function [ data, time, dt ] = load_2D_data( filename, variablename )
%LOAD_2D_DATA load a 2D variable stored in a hdf5 result file
    time     = h5read(filename,'/data/var2d/time');
    dt    = h5readatt(filename,'/data/input/basic','dt');
    cstart= h5readatt(filename,'/data/input/basic','start_iframe3d'); 
    
    % Find array size by loading the first output
    tmp   = h5read(filename,['/data/var2d/',variablename,'/', num2str(cstart+1,'%06d')]);
    try % check if it is complex or real
        sz  = size(tmp.real);
        cmpx = 1;
    catch
        sz  = size(tmp);
        cmpx = 0;
    end

    data     = zeros(sz);
    
    for it = 1:numel(time)
        tmp         = h5read(filename,['/data/var2d/',variablename,'/', num2str(cstart+it,'%06d')]);
        if cmpx
            data(:,:,it) = tmp.real + 1i * tmp.imaginary;
        else
            data(:,:,it) = tmp;
        end
    end

end

