function [TOPLOT] = poloidal_plot(DATA,OPTIONS)
sq = @(x) squeeze(x);
%% Geometry
r0 = DATA.a;
q0 = DATA.Q0;
Cy = r0/q0;
Ly = max(DATA.y) - min(DATA.y);
n0 = Cy*2*pi/Ly;
z_ = DATA.z;
%% grid sizes
nkx = DATA.Nkx; nky = DATA.Nky; nz = DATA.Nz;
nx  = DATA.Nx;  ny  = DATA.Ny;
%% Time and frames
FRAMES = zeros(size(OPTIONS.TIME));
for i = 1:numel(OPTIONS.TIME)
    [~,FRAMES(i)] =min(abs(OPTIONS.TIME(i)-DATA.Ts3D));
end
FRAMES = unique(FRAMES);

%% Field to plot
switch OPTIONS.NAME
    case '\phi'
        FIELD = DATA.PHI;
        NAME  = 'phi';
    case 'n_i'
        FIELD = DATA.DENS_I;
        NAME  = 'ni';
end

%% Build poloidal plot (tor angle = 0)
FRZ = zeros(DATA.Nx,DATA.Nz, numel(FRAMES));
w_ = zeros(nky,nz);
for iz = 1:nz
    for iky = 1:ny
        w_(iky,iz) = exp(1i*iky*n0*(q0*z_(iz)));
    end
end


for it = 1:numel(FRAMES)
    field_ = squeeze(FIELD(:,:,:,FRAMES(it)));
    spectrumKxKyZ = zeros(nkx,nky,nz);
    spectrumKxKyZ(1:nkx,:,:) = field_;
    spectrumKxKyZ((nkx+1):nx,1,:) = conj(field_(nkx:-1:2,1,:));
    spectrumKxKyZ((nkx+1):nx,2:ny,:) = conj(field_(nkx:-1:2,ny:-1:2,:));
    
    Fxkyz = ny*ifft(spectrumKxKyZ,[],1);
    for ix = 1:DATA.Nx
        FRZ(ix,:,it) = sq(sum(sq(Fxkyz(ix,:,:)).*w_(:,:),1));
    end
end
FRZ = real(FRZ);
%grid
[Y,X] = meshgrid(DATA.z,DATA.x);
X__ = (X+DATA.a).*cos(Y);
Y__ = (X+DATA.a).*sin(Y);
X = X__;
Y = Y__;
%% output

TOPLOT.FIELD     = FRZ;
TOPLOT.FRAMES    = FRAMES;
TOPLOT.X         = (X);
TOPLOT.Y         = (Y);
TOPLOT.FIELDNAME = NAME;
TOPLOT.XNAME     = '$R$';
TOPLOT.YNAME     = '$Z$';
TOPLOT.FILENAME  = [NAME,'_poloidal_phieq0'];
TOPLOT.DIMENSIONS= [100, 100, 500, 500];
TOPLOT.ASPECT    = [1 1 1];
TOPLOT.INTERP    = OPTIONS.INTERP;

end