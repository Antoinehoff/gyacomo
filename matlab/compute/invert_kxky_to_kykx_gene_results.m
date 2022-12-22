function [ data ] = invert_kxky_to_kykx_gene_results( data )
%to go from nkx nky Gene representation to nky nkx HeLaZ one

rotate = @(x) permute(x,[2 1 3 4]);

data.PHI    = rotate(data.PHI);
data.DENS_I = rotate(data.DENS_I);
data.TPER_I = rotate(data.TPER_I);
data.TPAR_I = rotate(data.TPAR_I);
data.TEMP_I = rotate(data.TEMP_I);

data.Ny = data.Nky*2-1;
data.Nx = data.Nkx;

if numel(data.kx)>1
    dkx = data.kx(2); 
else
    dkx = 1;
end

dky = data.ky(2);
Lx = 2*pi/dkx;   Ly = 2*pi/dky;
x = linspace(-Lx/2,Lx/2,data.Nx+1); x = x(1:end-1);
y = linspace(-Ly/2,Ly/2,data.Ny+1); y = y(1:end-1);
data.x = x; data.y = y; data.Lx = Lx; data.Ly = Ly;
end

