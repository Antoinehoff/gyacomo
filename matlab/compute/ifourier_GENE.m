function [ field_r ] = ifourier_GENE( field_c )
%IFOURIER_FIELD_GENE method of ifourier used in GENE post processing
%   This fourier transform back from the half complex plane to a full real
%   space, in 2 or 3D (dim). Compied from computeFXYZ of GENE for
%   comparison purpose.

%% Original
[nky,nx,nz]=size(field_c);
%extend to whole ky by imposing reality condition
ny=2*nky-1;

if ny~=1
    %note, we need one extra point which we set to zero for the ifft 
    spectrumKyKxZ=zeros(ny,nx,nz);
    spectrumKyKxZ(1:nky,:,:)=field_c(:,:,:);
    spectrumKyKxZ((nky+1):(ny),1,:)=conj(field_c(nky:-1:2,1,:));
    spectrumKyKxZ((nky+1):(ny),2:nx,:)=conj(field_c(nky:-1:2,nx:-1:2,:));
else
    %pad with zeros to interpolate on fine scale
    ny=20;
    spectrumKyKxZ=zeros(nx,ny,nz);
    spectrumKyKxZ(:,2,:)=field_c(:,:,:);
end   

%inverse fft, symmetric as we are using real data
spectrumXKyZ=nx*ifft(spectrumKyKxZ,[],1);
field_r=ny*ifft(spectrumXKyZ,[],2,'symmetric');
clear spectrumKxKyZ 
 
end

