function [ DATA ] = load_gene_data( folder )
%to load gene data as for HeLaZ results
namelist      = read_namelist([folder,'parameters']);
DATA.namelist = namelist;
DATA.folder   = folder;
%% Grid
coofile = 'coord.dat.h5';
DATA.grids.vp  = h5read([folder,coofile],'/coord/vp');
DATA.grids.Nvp = numel(DATA.grids.vp);

DATA.grids.mu  = h5read([folder,coofile],'/coord/mu');
DATA.grids.Nmu = numel(DATA.grids.mu);

DATA.grids.kx  = h5read([folder,coofile],'/coord/kx');
DATA.grids.Nkx = numel(DATA.grids.kx);
DATA.grids.Nx  = DATA.grids.Nkx;

DATA.grids.ky  = h5read([folder,coofile],'/coord/ky');
DATA.grids.Nky = numel(DATA.grids.ky);
DATA.grids.Ny  = DATA.grids.Nky*2-1;

DATA.grids.z   = h5read([folder,coofile],'/coord/z');
DATA.grids.Nz  = numel(DATA.grids.z);

DATA.params_string = [num2str(DATA.grids.Nkx),'x',num2str(DATA.grids.Nky),'x',num2str(DATA.grids.Nz),...
                    'x',num2str(DATA.grids.Nvp),'x',num2str(DATA.grids.Nmu)];
if numel(DATA.grids.kx)>1
    dkx = DATA.grids.kx(2); 
else
    dkx = 1;
end
dky = DATA.grids.ky(2);
Lx = 2*pi/dkx;   Ly = 2*pi/dky;
x = linspace(-Lx/2,Lx/2,DATA.grids.Nx+1); x = x(1:end-1);
y = linspace(-Ly/2,Ly/2,DATA.grids.Ny+1); y = y(1:end-1);
DATA.grids.x = x; DATA.grids.y = y; DATA.grids.Lx = Lx; DATA.grids.Ly = Ly;
%% Transport
nrgfile           = 'nrg.dat.h5';
% nrgfile           = 'nrg_1.h5';
DATA.Ts0D      = h5read([folder,nrgfile],'/nrgions/time');
DATA.PGAMMA_RI = h5read([folder,nrgfile],'/nrgions/Gamma_es');
DATA.HFLUX_X   = h5read([folder,nrgfile],'/nrgions/Q_es');

%% fields and moments
phifile   = 'field.dat.h5';
% phifile   = 'field_1.h5';
DATA.Ts3D = h5read([folder,phifile],'/field/time');
DATA.DENS_I = zeros(DATA.grids.Nkx,DATA.grids.Nky,DATA.grids.Nz,numel(DATA.Ts3D));
DATA.TPER_I = zeros(DATA.grids.Nkx,DATA.grids.Nky,DATA.grids.Nz,numel(DATA.Ts3D));
DATA.TPAR_I = zeros(DATA.grids.Nkx,DATA.grids.Nky,DATA.grids.Nz,numel(DATA.Ts3D));
DATA.PHI    = zeros(DATA.grids.Nkx,DATA.grids.Nky,DATA.grids.Nz,numel(DATA.Ts3D));

momfile = 'mom_ions.dat.h5';
% momfile = 'mom_ions_1.h5';
for jt = 1:numel(DATA.Ts3D)
    t = DATA.Ts3D(jt);
    [~, it] = min(abs(DATA.Ts3D-t));
try
    tmp = h5read([folder,momfile],['/mom_ions/dens/',sprintf('%10.10d',it-1)]);
    DATA.DENS_I(:,:,:,it) = tmp.real + 1i*tmp.imaginary;
catch
    DATA.DENS_I(:,:,:,it) = 0;
end
try
    tmp = h5read([folder,momfile],['/mom_ions/T_par/',sprintf('%10.10d',it-1)]);
    DATA.TPAR_I(:,:,:,it) = tmp.real + 1i*tmp.imaginary;
catch
    DATA.TPAR_I(:,:,:,it) = 0;
end
try
    tmp = h5read([folder,momfile],['/mom_ions/T_perp/',sprintf('%10.10d',it-1)]);
    DATA.TPER_I(:,:,:,it) = tmp.real + 1i*tmp.imaginary;
catch
    DATA.TPER_I(:,:,:,it) = 0;
end
    tmp = h5read([folder,phifile],['/field/phi/',sprintf('%10.10d',it-1)]);
    DATA.PHI(:,:,:,it) = tmp.real + 1i*tmp.imaginary;

end

DATA.TEMP_I = (DATA.TPAR_I + 2*DATA.TPER_I)/3.0 - DATA.DENS_I;
DATA.scale = 1;
DATA.PARAMS = ['GENE'];
DATA.param_title = 'GENE';
DATA.localdir = folder;
%% Geometry
CMD = ['tail -n ',num2str(namelist.box.nz0),' ',folder,namelist.geometry.magn_geometry{1},'.dat > tmp.txt'];
system(CMD);
DATA.geo_arrays = load('tmp.txt');
system('rm tmp.txt');
end

