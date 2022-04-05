function [ data ] = rotate_c_plane_nxnky_to_nkxny( data )
%to go from nx nky Gene representation to nkx ny HeLaZ one

kx_full  = data.kx;

ky_full  = [data.ky; -data.ky(end-1:-1:2)];

dens = data.DENS_I; temp = data.TEMP_I; phi = data.PHI;

dims = size(dens); nz = dims(3); nt = dims(4);

nkx = numel(data.kx); 
nky = numel(data.ky); 
nx  = nkx;
ny  = 2*numel(data.ky)-1;

%note, we need one extra point which we set to zero for the ifft 
dens_full                        = zeros(nx,ny,nz,nt);
dens_full(:,1:nky,:,:)           = dens(:,:,:,:);
dens_full(1,(nky+1):(ny),:,:)    = conj(dens(1,nky:-1:2,:,:));       
dens_full(2:nx,(nky+1):(ny),:,:) = conj(dens(nx:-1:2,nky:-1:2,:,:)); 
dens_full = cat(2,dens(1:data.Nkx,:,:,:),conj(dens(1:data.Nkx,end-1:-1:2,:,:)));

temp_full                        = zeros(nx,ny,nz,nt);
temp_full(:,1:nky,:,:)           = temp(:,:,:,:);
temp_full(1,(nky+1):(ny),:,:)    = conj(temp(1,nky:-1:2,:,:));       
temp_full(2:nx,(nky+1):(ny),:,:) = conj(temp(nx:-1:2,nky:-1:2,:,:)); 
temp_full = cat(2,temp(1:data.Nkx,:,:,:),conj(temp(1:data.Nkx,end-1:-1:2,:,:)));

phi_full                        = zeros(nx,ny,nz,nt);
phi_full(:,1:nky,:,:)           = phi(:,:,:,:);
phi_full(1,(nky+1):(ny),:,:)    = conj(phi(1,nky:-1:2,:,:));       
phi_full(2:nx,(nky+1):(ny),:,:) = conj(phi(nx:-1:2,nky:-1:2,:,:));
phi_full  = cat(2, phi(1:data.Nkx,:,:,:),conj( phi(1:data.Nkx,end-1:-1:2,:,:)));

data.DENS_I = dens_full(1:data.Nkx/2+1,:,:,:);
data.TEMP_I = temp_full(1:data.Nkx/2+1,:,:,:);
data.PHI    = phi_full (1:data.Nkx/2+1,:,:,:);
data.kx     = kx_full(1:data.Nkx/2+1);
data.ky     = ky_full;
data.Nkx    = numel(data.kx);
data.Nky    = numel(data.ky);
data.Nx     = nx+1;
data.Ny     = ny-1;

dkx = data.kx(2); dky = data.ky(2);
Lx = 2*pi/dkx;   Ly = 2*pi/dky;
x = linspace(-Lx/2,Lx/2,data.Nx+1); x = x(1:end-1);
y = linspace(-Ly/2,Ly/2,data.Ny+1); y = y(1:end-1);
data.x = x; data.y = y; data.Lx = Lx; data.Ly = Ly;
end

