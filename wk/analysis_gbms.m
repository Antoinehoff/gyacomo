% cd /home/ahoffman/Documents/gbms/scan_test
% addpath(genpath('/home/ahoffman/Documents/gbms/matlab_scripts'));
% res = gbms_get_scandir('/home/ahoffman/Documents/gbms/scan_test/scan_test/');
% figure; plot(res.paramscan,res.growth_rate)

%%
% resdir = '/home/ahoffman/Documents/gbms/benchmark_HeLaZ/shearless_linear_cyclone/';
% resdir = '/home/ahoffman/Documents/gbms/benchmark_HeLaZ/new_RH_test/';
% resdir = '/home/ahoffman/Documents/gbms/benchmark_HeLaZ/RH_test_kine/';
% resdir = '/home/ahoffman/Documents/gbms/benchmark_HeLaZ/KBM/';
% resdir = '/home/ahoffman/Documents/gbms/benchmark_HeLaZ/TEM/';
% resdir = '/home/ahoffman/Documents/gbms/benchmark_HeLaZ/ITG/';
% resdir = '/home/ahoffman/Documents/gbms/benchmark_HeLaZ/linear_cyclone/';
resdir = '/home/ahoffman/Documents/gbms/results/MTM/';

% system(['cd ',resdir,';','./gbms < parameters.in; cd /home/ahoffman/HeLaZ/wk']);
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
gbms_dat.BETA    = h5readatt(outfile,'/data/input','betae');
gbms_dat.SHEAR    = h5readatt(outfile,'/data/input','magnetic shear');
gbms_dat.PHI = zeros(gbms_dat.Ny,gbms_dat.Nx,gbms_dat.Nz,gbms_dat.Nt);
gbms_dat.PSI = zeros(gbms_dat.Ny,gbms_dat.Nx,gbms_dat.Nz,gbms_dat.Nt);
gbms_dat.param_title = 'GBMS';
for it = 1:gbms_dat.Nt
    
    tmp = h5read(outfile,['/data/var2d/phi/',sprintf('%.6d',it-1)]);
    gbms_dat.PHI(:,:,:,it) = permute(tmp.real + 1i * tmp.imaginary,[2 1 3]);
    if gbms_dat.BETA > 0
        tmp = h5read(outfile,['/data/var2d/psi/',sprintf('%.6d',it-1)]);
        gbms_dat.PSI(:,:,:,it) = permute(tmp.real + 1i * tmp.imaginary,[2 1 3]);
    end

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
trange = [0.5 1]*gbms_dat.Ts3D(end);
nplots = 3;
lg = compute_fluxtube_growth_rate(gbms_dat,trange,nplots);
[gmax,     kmax] = max(lg.g_ky(:,end));
[gmaxok, kmaxok] = max(lg.g_ky(:,end)./lg.ky);
msg = sprintf('gmax = %2.2f, kmax = %2.2f',gmax,lg.ky(kmax)); disp(msg);
msg = sprintf('gmax/k = %2.2f, kmax/k = %2.2f',gmaxok,lg.ky(kmaxok)); disp(msg);
end

if 1
%% Ballooning plot
% options.time_2_plot = data.Ts3D(end);
% options.kymodes     = [0.5];
% options.normalized  = 1;
% options.sheared     = 0;
% options.field       = 'phi';
% fig = plot_ballooning(gbms_dat,options);
plot_gbms_ballooning(outfile);

% plot(b_angle,phib_real); hold on;

end

if 0
%% RH TEST
ikx = 1; iky = 1; t0 = 0; t1 = gbms_dat.Ts3D(end);
[~, it0] = min(abs(t0-gbms_dat.Ts3D));[~, it1] = min(abs(t1-gbms_dat.Ts3D));
plt = @(x) squeeze(mean(real(x(iky,ikx,:,it0:it1)),3));%./squeeze(mean(real(x(iky,ikx,:,it0)),3));
figure
plot(gbms_dat.Ts3D(it0:it1), plt(gbms_dat.PHI),'k');
xlabel('$t$'); ylabel('$\phi_z(t)/\phi_z(0)$')
title(sprintf('$k_x=$%2.2f, $k_y=$%2.2f',gbms_dat.kx(ikx),gbms_dat.ky(iky)))
end