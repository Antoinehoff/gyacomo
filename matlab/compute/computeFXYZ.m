function output=computeFXYZ(gene_data,varargin)
%Function for antitrasforming a Fourier field to real space
% 
%  INPUTS:
%      gene_data  -> 3D field in Fourier
%      x_local    -> Fourier in x (1st dim.)
%
%    [output]=COMPUTEFXYZ(gene_data)

%%
if nargin==1
    x_local=1;
else
    x_local=varargin{1};
end

if x_local
    [nx,nky,nz]=size(gene_data);
    %extend to whole ky by imposing reality condition
    ny=2*nky-1;
    
    if ny~=1
        %note, we need one extra point which we set to zero for the ifft 
        spectrumKxKyZ=zeros(nx,ny,nz);
        spectrumKxKyZ(:,1:nky,:)=gene_data(:,:,:);
        spectrumKxKyZ(1,(nky+1):(ny),:)=conj(gene_data(1,nky:-1:2,:));
        spectrumKxKyZ(2:nx,(nky+1):(ny),:)=conj(gene_data(nx:-1:2,nky:-1:2,:));
    else
        %pad with zeros to interpolate on fine scale
        ny=20;
        spectrumKxKyZ=zeros(nx,ny,nz);
        spectrumKxKyZ(:,2,:)=gene_data(:,:,:);
    end   
    %inverse fft, symmetric as we are using real data
    spectrumXKyZ=nx*ifft(spectrumKxKyZ,[],1);
    output=ny*ifft(spectrumXKyZ,[],2,'symmetric');
    clear  spectrumKxKyZ 
 
else  
    
    [nx,nky,nz]=size(gene_data);
    %   extend to whole ky by imposingreality condition
    ny=2*nky-1;
    if ny~=1
        spectrumXKyZ=zeros(nx,ny,nz);
        spectrumXKyZ(:,1:nky,:)=gene_data(:,:,:);
        spectrumXKyZ(1,(nky+1):(ny),:)=conj(gene_data(1,nky:-1:2,:));
        spectrumXKyZ(2:nx,(nky+1):(ny),:)=conj(gene_data(2:1:nx,nky:-1:2,:));   
    else
        %pad with zeros to interpolate on fine scale
        ny=50;
        spectrumXKyZ=zeros(nx,ny,nz);
        spectrumXKyZ(:,2,:)=gene_data(:,:,:);

    end
    output=ny*ifft(spectrumXKyZ,[],2,'symmetric');
    clear  spectrumK=XKyZ 
end
