function [ data, time, dt ] = load_pjz_data( filename,variablename,specie)
    time     = h5read(filename,'/data/var3d/time');
    p        = h5read(filename,['/data/grid/coordp_',specie]);
    j        = h5read(filename,['/data/grid/coordj_',specie]);
    z        = h5read(filename,'/data/grid/coordz');
    dt    = h5readatt(filename,'/data/input','dt');
    cstart= h5readatt(filename,'/data/input','start_iframe3d'); 
    
    data     = zeros(numel(p),numel(j),numel(z),numel(time));
    for it = 1:numel(time)
        tmp         = h5read(filename,['/data/var3d/',variablename,'/', num2str(cstart+it,'%06d')]);
        data(:,:,:,it) = tmp;
    end

end

