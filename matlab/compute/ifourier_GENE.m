function [ field_r ] = ifourier_GENE( field_c, dims )
%IFOURIER_FIELD_GENE method of ifourier used in GENE post processing
%   This fourier transform back from the half complex plane to a full real
%   space, in 2 or 3D (dim). Compied from computeFXYZ of GENE for
%   comparison purpose.
sz = size(field_c);
nkx = sz(1);
field_r = zeros(dims);

if numel(dims) == 2 %2D
    nx = dims(1); ny = dims(2);
    %note, we need one extra point which we set to zero for the ifft 
    spectrumKxKyZ                  = zeros(nx,ny);
    spectrumKxKyZ(1:nkx,:)         = field_c;
    spectrumKxKyZ((nkx):(nx),1)    = conj(field_c(nkx:-1:2));       
    spectrumKxKyZ((nkx):(nx),2:ny) = conj(field_c(nkx:-1:2,ny:-1:2));
    
%     spectrumKxYZ = ny*ifft(spectrumKxKyZ,[],2);
%     field_r      = nx*ifft(spectrumKxYZ,[],1,'symmetric');
    
    spectrumKxYZ = ifft(spectrumKxKyZ,[],2);
    field_r      = ifft(spectrumKxYZ,[],1,'symmetric');
end

end

