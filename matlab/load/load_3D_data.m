function [ data, time, dt ] = load_3D_data( filename, variablename )
%LOAD_2D_DATA load a 2D variable stored in a hdf5 result file from HeLaZ
    time     = h5read(filename,'/data/var3d/time');
    kx       = h5read(filename,'/data/grid/coordkx');
    ky       = h5read(filename,'/data/grid/coordky');
    z        = h5read(filename,'/data/grid/coordz');
    dt    = h5readatt(filename,'/data/input','dt');
    cstart= h5readatt(filename,'/data/input','start_iframe3d'); 
    
    data     = zeros(numel(ky),numel(kx),numel(z),numel(time));
    [KX,KY]  = meshgrid(kx,ky);
    
    switch variablename
    case 'vEx'
         for it = 1:numel(time)
            tmp         = h5read(filename,['/data/var3d/phi/', num2str(cstart+it,'%06d')]);
            data(:,:,:,it) = 1i.*KY.*(tmp.real + 1i * tmp.imaginary);
         end
    case 'vEy'
         for it = 1:numel(time)
            tmp         = h5read(filename,['/data/var3d/phi/', num2str(cstart+it,'%06d')]);
            data(:,:,:,it) = -1i.*KX.*(tmp.real + 1i * tmp.imaginary);
        end           
    otherwise
        for it = 1:numel(time)
            tmp         = h5read(filename,['/data/var3d/',variablename,'/', num2str(cstart+it,'%06d')]);
            data(:,:,:,it) = tmp.real + 1i * tmp.imaginary;
        end
    end

end

