function [ field_r ] = ifourier_GENE( field_c )
%IFOURIER_FIELD_GENE method of ifourier used in GENE post processing
%   This fourier transform back from the half complex plane to a full real
%   space, in 2 or 3D (dim). Compied from computeFXYZ of GENE for
%   comparison purpose.

%% Original
[nx,nky,nz]=size(field_c);
%extend to whole ky by imposing reality condition
ny=2*nky-1;

if ny~=1
    %note, we need one extra point which we set to zero for the ifft 
    spectrumKxKyZ=zeros(nx,ny,nz);
    spectrumKxKyZ(:,1:nky,:)=field_c(:,:,:);
    spectrumKxKyZ(1,(nky+1):(ny),:)=conj(field_c(1,nky:-1:2,:));
    spectrumKxKyZ(2:nx,(nky+1):(ny),:)=conj(field_c(nx:-1:2,nky:-1:2,:));
else
    %pad with zeros to interpolate on fine scale
    ny=20;
    spectrumKxKyZ=zeros(nx,ny,nz);
    spectrumKxKyZ(:,2,:)=field_c(:,:,:);
end   

%inverse fft, symmetric as we are using real data
spectrumXKyZ=nx*ifft(spectrumKxKyZ,[],1);
field_r=ny*ifft(spectrumXKyZ,[],2,'symmetric');
clear spectrumKxKyZ 

%% Adapted for HeLaZ old representation
% [nkx,ny,nz]=size(field_c);
% %extend to whole ky by imposing reality condition
% nx=2*nkx-1;
% 
% %note, we need one extra point which we set to zero for the ifft 
% spectrumKxKyZ=zeros(nx,ny,nz);
% spectrumKxKyZ(1:nkx,:,:)=field_c(:,:,:);
% spectrumKxKyZ((nkx+1):(nx),1,:)=conj(field_c(nkx:-1:2,1,:));
% spectrumKxKyZ((nkx+1):(nx),2:ny,:)=conj(field_c(nkx:-1:2,ny:-1:2,:));
% 
% %inverse fft, symmetric as we are using real data
% spectrumXKyZ=nx*ifft(spectrumKxKyZ,[],1);
% field_r=ny*ifft(spectrumXKyZ,[],2,'symmetric');
% clear spectrumKxKyZ 
 
end

