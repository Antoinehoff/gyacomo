function [ data ] = invert_kxky_to_kykx_gene_results( data )
%to go from nkx nky Gene representation to nky nkx HeLaZ one

rotate = @(x) permute(x,[2 1 3 4]);

data.PHI    = rotate(data.PHI);
data.DENS_I = rotate(data.DENS_I);
data.TPER_I = rotate(data.TPER_I);
data.TPAR_I = rotate(data.TPAR_I);
data.TEMP_I = rotate(data.TEMP_I);

data.grids.Ny = data.grids.Nky*2-1;
data.grids.Nx = data.grids.Nkx;

if numel(data.grids.kx)>1
    dkx = data.grids.kx(2); 
else
    dkx = 1;
end
if data.grids.Nky > 1
    dky = data.grids.ky(2);
else
    dky = data.grids.ky(1);
end
Lx = 2*pi/dkx;   Ly = 2*pi/dky;
x = linspace(-Lx/2,Lx/2,data.grids.Nx+1); x = x(1:end-1);
y = linspace(-Ly/2,Ly/2,data.grids.Ny+1); y = y(1:end-1);
data.grids.x = x; data.grids.y = y; data.grids.Lx = Lx; data.grids.Ly = Ly;
end

