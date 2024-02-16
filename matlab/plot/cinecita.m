function [ ] = cinecita(DATADIR,J0,J1,fieldname,plan,save)

data   = {};
[data] = compile_results_low_mem(data,DATADIR,J0,J1);
%% Select the field
SKIP_COMP = 0; % Turned on only for 2D plots
OPE_      = 1; % Operation on data
switch fieldname
    case 'phi' %ES pot
        [FIELD,TIME] = compile_results_3D(DATADIR,J0,J1,'phi');
        ltxname = '\phi';
    case 'phi_obmp' %ES pot
        [FIELD,TIME] = compile_results_2D(DATADIR,J0,J1,'phi_obmp');
        ltxname = '\phi_{z=0}';
        SKIP_COMP = 1;
    case '\psi' %EM pot
        [FIELD,TIME] = compile_results_3D(DATADIR,J0,J1,'psi');
        ltxname = '\psi';
    case 'phi_nz' % non-zonal ES pot
        [FIELD,TIME] = compile_results_3D(DATADIR,J0,J1,'phi');
        ltxname = '\phi^{NZ}';
        OPE_ = (KY~=0);  
    case 'vE_y' % y-comp of ExB velocity
        [FIELD,TIME] = compile_results_3D(DATADIR,J0,J1,'phi');
        ltxname = 'v_{Ey}';
        OPE_ = -1i*KX;  
   case 'vE_x' % x-comp of ExB velocity
        [FIELD,TIME] = compile_results_3D(DATADIR,J0,J1,'phi');
        ltxname = 'v_{Ex}';
        OPE_ = -1i*KY;  
   case 'sE_y' % y-comp of ExB shear
        [FIELD,TIME] = compile_results_3D(DATADIR,J0,J1,'phi');
        ltxname = 's_{Ey}';
        OPE_ = KX.^2; 
   case 'sE_x' % x-comp of ExB shear
        [FIELD,TIME] = compile_results_3D(DATADIR,J0,J1,'phi');
        ltxname = 's_{Ex}';
        OPE_ = KY.^2; 
   case 'wz' % ES pot vorticity
        [FIELD,TIME] = compile_results_3D(DATADIR,J0,J1,'phi');
        ltxname = '\omega_z';
        OPE_ = -(KX.^2+KY.^2);        
   otherwise
        disp('Fieldname not recognized');
        return
end
%% Setup  grids
% [KX, KY] = meshgrid(DATA.grids.kx, DATA.grids.ky);
Nx = data.grids.Nx; Ny = data.grids.Ny; Nz = data.grids.Nz; Nt = numel(TIME);
% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
options.INTERP    = 0;
options.POLARPLOT = 0;
options.BWR       = 1; % bluewhitered plot or gray
options.CLIMAUTO  = 1; % adjust the colormap auto
options.NAME      = ltxname;
options.PLAN      = plan;
options.COMP      = 1;
options.TIME      =  data.Time(1:1:end);
data.EPS          = 0.1;
data.a = data.EPS * 2000;
options.RESOLUTION = 256;
create_film(data,options,'.gif')

if ~save
    system(['rm ',FILENAME]);
end

end