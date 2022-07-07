function [ pola ] = compute_polarisation( data )
%compute_polarisation computes the polarisation term (1-Gamma0)phi with
%ordering up to Jmax
%   Detailed explanation goes here
PHI = data.PHI;
T   = data.Ts3D;
Ji  = data.Jmaxi;
pola = zeros(size(PHI));

KN_ = zeros(data.Nky,data.Nkx,data.Nz);

[KX,KY] = meshgrid(data.kx,data.ky);
KP  = sqrt(KX.^2+KY.^2);

GAMMA2_ = 0.*KN_;
   
for in = 1:Ji+1
    for iz = 1:data.Nz
        GAMMA2_(:,:,iz) = GAMMA2_(:,:,iz) + kernel(in-1,KP*sqrt(2)).^2;
    end
end


for it = 1:numel(T)
    
    pola(:,:,:,it) = (1-GAMMA2_).*PHI(:,:,:,it);
    
end


end

