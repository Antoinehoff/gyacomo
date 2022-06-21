% cd /home/ahoffman/Documents/gbms/scan_test
% addpath(genpath('/home/ahoffman/Documents/gbms/matlab_scripts'));
% res = gbms_get_scandir('/home/ahoffman/Documents/gbms/scan_test/scan_test/');
% figure; plot(res.paramscan,res.growth_rate)

%%
% resdir = '/home/ahoffman/Documents/gbms/benchmark_HeLaZ/shearless_linear_cyclone/';
resdir = '/home/ahoffman/Documents/gbms/benchmark_HeLaZ/RH_test/';
% resdir = '/home/ahoffman/Documents/gbms/benchmark_HeLaZ/linear_cyclone/';
% resdir = '/home/ahoffman/molix/';
outfile = [resdir,'field.dat.h5'];

gbms_dat.Ts3D    = h5read(outfile,'/data/var2d/time');
gbms_dat.Nt      = unique(numel(gbms_dat.Ts3D));
gbms_dat.kx      = unique(h5read(outfile,'/data/var2d/phi/coordkx'));
gbms_dat.ky      = unique(h5read(outfile,'/data/var2d/phi/coordky'));
gbms_dat.z       = unique(h5read(outfile,'/data/var2d/phi/coordz'));
gbms_dat.Nx = numel(gbms_dat.kx); gbms_dat.Nkx = numel(gbms_dat.kx); 
gbms_dat.Ny = numel(gbms_dat.ky); gbms_dat.Nky = numel(gbms_dat.ky); 
gbms_dat.Nz = numel(gbms_dat.z);

dky = min(gbms_dat.ky(gbms_dat.ky>0)); Ly =0;% 2*pi/dky;
gbms_dat.y  = linspace(-Ly/2,Ly/2,gbms_dat.Ny+1); gbms_dat.y = gbms_dat.y(1:end-1);
gbms_dat.x  = 0;
gbms_dat.PHI = zeros(gbms_dat.Ny,gbms_dat.Nx,gbms_dat.Nz,gbms_dat.Nt);
gbms_dat.param_title = 'GBMS';
for it = 1:gbms_dat.Nt
    
    tmp = h5read(outfile,['/data/var2d/phi/',sprintf('%.6d',it-1)]);
    gbms_dat.PHI(:,:,:,it) = permute(tmp.real + 1i * tmp.imaginary,[2 1 3]);
    
end

gbms_dat.localdir = resdir;

%%
if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
options.INTERP    = 1;
options.POLARPLOT = 0;
options.NAME      = '\phi';
options.PLAN      = 'yz';
options.COMP      = 1;
options.TIME      = 0:200;
gbms_dat.EPS          = 0.1
gbms_dat.a = gbms_dat.EPS * 2000;
create_film(gbms_dat,options,'.gif')
end

if 0
%% 2D snapshots
% Options
options.INTERP    = 0;  
options.POLARPLOT = 0;
options.AXISEQUAL = 1;
options.NAME      = '\phi';
options.PLAN      = 'yz';
options.COMP      = 1;
options.TIME      = 100;
gbms_dat.EPS = 1e-3;
gbms_dat.a = gbms_dat.EPS * 2000;
fig = photomaton(gbms_dat,options);
save_figure(gbms_dat,fig)
end

if 0
%% linear growth rate for 3D fluxtube
trange = [10 200];
nplots = 1;
lg = compute_fluxtube_growth_rate(gbms_dat,trange,nplots);
end

if 0
%% Ballooning plot
options.time_2_plot = data.Ts3D(end);
options.kymodes     = [0.5];
options.normalized  = 1;
options.sheared     = 0;
options.field       = 'phi';
fig = plot_ballooning(gbms_dat,options);
end

if 1
%% RH TEST
ikx = 1;
plt = @(x) squeeze(mean(real(x(1,ikx,:,:)),3))./squeeze(mean(real(x(1,ikx,:,1)),3));
figure
plot(gbms_dat.Ts3D, plt(gbms_dat.PHI));
xlabel('$t$'); ylabel('$\phi_z(t)/\phi_z(0)$')
title(sprintf('$k_x=$%2.2f, $k_y=0.00$',gbms_dat.kx(ikx)))
end