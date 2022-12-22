function [ DATA ] = load_gene_data( folder )
%to load gene data as for HeLaZ results
namelist      = read_namelist([folder,'parameters.dat']);
DATA.namelist = namelist;
%% Grid
coofile = 'coord.dat.h5';
DATA.vp  = h5read([folder,coofile],'/coord/vp');
DATA.Nvp = numel(DATA.vp);

DATA.mu  = h5read([folder,coofile],'/coord/mu');
DATA.Nmu = numel(DATA.mu);

DATA.kx  = h5read([folder,coofile],'/coord/kx');
DATA.Nkx = numel(DATA.kx);
DATA.Nx  = DATA.Nkx;

DATA.ky  = h5read([folder,coofile],'/coord/ky');
DATA.Nky = numel(DATA.ky);
DATA.Ny  = DATA.Nky*2-1;

DATA.z   = h5read([folder,coofile],'/coord/z');
DATA.Nz  = numel(DATA.z);

if numel(DATA.kx)>1
    dkx = DATA.kx(2); 
else
    dkx = 1;
end
dky = DATA.ky(2);
Lx = 2*pi/dkx;   Ly = 2*pi/dky;
x = linspace(-Lx/2,Lx/2,DATA.Nx+1); x = x(1:end-1);
y = linspace(-Ly/2,Ly/2,DATA.Ny+1); y = y(1:end-1);
DATA.x = x; DATA.y = y; DATA.Lx = Lx; DATA.Ly = Ly;
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
DATA.DENS_I = zeros(DATA.Nkx,DATA.Nky,DATA.Nz,numel(DATA.Ts3D));
DATA.TPER_I = zeros(DATA.Nkx,DATA.Nky,DATA.Nz,numel(DATA.Ts3D));
DATA.TPAR_I = zeros(DATA.Nkx,DATA.Nky,DATA.Nz,numel(DATA.Ts3D));
DATA.PHI    = zeros(DATA.Nkx,DATA.Nky,DATA.Nz,numel(DATA.Ts3D));

momfile = 'mom_ions.dat.h5';
% momfile = 'mom_ions_1.h5';
for jt = 1:numel(DATA.Ts3D)
    t = DATA.Ts3D(jt);
    [~, it] = min(abs(DATA.Ts3D-t));
% 
%     tmp = h5read([folder,momfile],['/mom_ions/dens/',sprintf('%10.10d',it-1)]);
%     DATA.DENS_I(:,:,:,it) = tmp.real + 1i*tmp.imaginary;
% % 
%     tmp = h5read([folder,momfile],['/mom_ions/T_par/',sprintf('%10.10d',it-1)]);
%     DATA.TPAR_I(:,:,:,it) = tmp.real + 1i*tmp.imaginary;
% %  
%     tmp = h5read([folder,momfile],['/mom_ions/T_perp/',sprintf('%10.10d',it-1)]);
%     DATA.TPER_I(:,:,:,it) = tmp.real + 1i*tmp.imaginary;
%     
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

