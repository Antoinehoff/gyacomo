% folder = '/misc/gene_results/NL_Zpinch_Kn_1.8_eta_0.25_nuSG_5e-2_mu_1e-2_KHR_2/';
% folder = '/misc/gene_results/NL_Zpinch_Kn_1.8_eta_0.25_nuSG_5e-3_gyroLES/';
% folder = '/misc/gene_results/HP_fig_2c_gyroLES/';
folder = '/misc/gene_results/NL_Zpinch_Kn_1.8_eta_0.25_nuSG_5e-2_mu_1e-2_SGDK/';

%% General info (domain etc.)
varargin{1} = 0; varargin{2} = 0;
general = loadGeneral_GENE(folder,-1,varargin);

nx = general.allParams.box.nx; Lx = general.allParams.box.Lx;
ny = general.allParams.box.ny; Ly = general.allParams.box.Ly;
nz = general.allParams.box.nz;
kx = general.coord.kx;
ky = general.coord.ky;
[KY,KX] = meshgrid(ky,kx);
x  = Lx/nx.*(-nx/2:nx/2-1);
y  = Ly/ny.*(-ny+1:ny-1);

%% Load dens_i
steps_range = 0:1:500; nt = numel(steps_range);
time_dens        = 1:nt;
dens_xt = zeros(nt,nx);
% figure
for i = 1:nt
    field_step = steps_range(i);
    
    [field_c, t] = loadMomentumCpx(folder,-1,field_step-1,'ions','dens','none',0,[nx ny nz],0,general.allParams.data_format);
    
    field_r = computeFXYZ(field_c,1);
    dens_xt(i,:) = mean(field_r(:,:,1),2);
    time_dens(i) = t;
end

if 1
%% Load phi
steps_range = 0:2:1000; nt = numel(steps_range);
time_phi        = 1:nt;
phi_xt = zeros(nt,nx);
% figure
for i = 1:nt
    field_step = steps_range(i);
    
    [field_c, t] = loadFieldPhiCpx(folder,-1,field_step-1,'',0,0,[nx ny nz],general.allParams.data_format);
    
    field_r = computeFXYZ(field_c,1);
    phi_xt(i,:) = mean(field_r(:,:,1),2);
    time_phi(i) = t;
end
end

if 0
%% Load v_ExB
steps_range = 0:2:1000; nt = numel(steps_range);
time_vExB        = 1:nt;
vExB_xt = zeros(nt,nx);
% figure
for i = 1:nt
    field_step = steps_range(i);
    
    [field_c, t] = loadFieldPhiCpx(folder,-1,field_step-1,'',0,0,[nx ny nz],general.allParams.data_format);
    
    field_r = computeFXYZ(-KX.*field_c,1);
    vExB_xt(i,:) = mean(field_r(:,:,1),2);
    time_vExB(i) = t;
end
end
if 1
%% compute Gx
steps_n = 0:1:500;  nt1 = numel(steps_n);
steps_v = 0:2:1000; nt2 = numel(steps_v);
time_vExB        = 1:nt1;
Gx = zeros(nt,1);
% figure
for i = 1:nt1
    field_step = steps_n(i);
    [field_c, t] = loadMomentumCpx(folder,-1,field_step-1,'ions','dens','none',0,[nx ny nz],0,general.allParams.data_format);
    dens = computeFXYZ(field_c,1);
    time_dens(i) = t;
    
    field_step = steps_v(i);
    [field_c, t] = loadFieldPhiCpx(folder,-1,field_step-1,'',0,0,[nx ny nz],general.allParams.data_format);
    vExB = computeFXYZ(-KY.*field_c,1);
    time_vExB(i) = t;
    
    Gx(i) = mean(mean(squeeze(dens(:,:,1).*vExB(:,:,1))));
end
end
%% Plot
figure; 
[X_TX,Y_TX] = meshgrid(time_dens,x);
% pclr=pcolor(X_TX',Y_TX',dens_xt);        
pclr=pcolor(X_TX',Y_TX',phi_xt);        
set(pclr, 'edgecolor','none'); 
colormap(bluewhitered(256))

%% 
% figure
% plot(time_dens,Gx)
