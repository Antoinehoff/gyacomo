addpath(genpath('../matlab')) % ... add
addpath(genpath('../matlab/plots')) % ... add
%% Directory of the simulation
if 1% Local results
outfile ='';
outfile ='';
outfile ='';
outfile ='';
outfile ='HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e-02_DGGK_mu_3e-03';
    BASIC.RESDIR      = ['../results/',outfile,'/'];
    BASIC.MISCDIR     = ['/misc/HeLaZ_outputs/results/',outfile,'/'];
    CMD = ['cp ', BASIC.RESDIR,'outputs* ',BASIC.MISCDIR]; disp(CMD);
    system(CMD);
else% Marconi results
outfile ='';
outfile ='';
outfile ='';
outfile ='';
outfile ='';
outfile ='';
% outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e-01_SGGK_mu_0e+00/out.txt';
outfile ='/marconi_scratch/userexternal/ahoffman/HeLaZ/results/HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e-01_SGGK_mu_3e-02/out.txt';
BASIC.RESDIR      = ['../',outfile(46:end-8),'/'];
BASIC.MISCDIR     = ['/misc/HeLaZ_outputs/',outfile(46:end-8),'/'];
end

%% Load the results
% Load outputs from jobnummin up to jobnummax
JOBNUMMIN = 00; JOBNUMMAX = 20; 
compile_results %Compile the results from first output found to JOBNUMMAX if existing

%% Post-processing
post_processing

%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default_plots_options
disp('Plots')
FMT = '.fig';

if 0
%% Time evolutions and growth rate
plot_time_evolution_and_gr
end

if 1
%% Space time diagramm (fig 11 Ivanov 2020)
TAVG = 1000; % Averaging time duration
plot_radial_transport_and_shear
end

if 0
%% Space time diagramms
cmax = 0.01 % max of the colorbar for transport
tstart = 0; tend = Ts3D(end); % time window
plot_space_time_diagrams
end

if 0
%% |phi_k|^2 spectra (Kobayashi 2015 fig 3)
% tstart = 0.8*Ts3D(end); tend = Ts3D(end); % Time window
tstart = 5000; tend = tstart+1000;
% Chose the field to plot
% FIELD = Ni00;   FNAME = 'Ni00'; FIELDLTX = 'N_i^{00}';
% FIELD = Ne00;   FNAME = 'Ne00'; FIELDLTX = 'N_e^{00}'
% FIELD = PHI; FNAME = 'PHI'; FIELDLTX = '\tilde\phi';
FIELD_ = fft2(Gamma_x); FIELD = FIELD_(1:76,:,:,:); FNAME = 'Gamma_x'; FIELDLTX = '\tilde\Gamma_x';
LOGSCALE = 1; TRENDS = 0;
plot_kperp_spectrum
end

if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
t0    =3000; iz = 1; ix = 1; iy = 1;
skip_ =1; DELAY = 1e-2*skip_;
[~, it03D] = min(abs(Ts3D-t0)); FRAMES_3D = it03D:skip_:numel(Ts3D);
[~, it05D] = min(abs(Ts5D-t0)); FRAMES_5D = it05D:skip_:numel(Ts5D);
INTERP = 0; T = Ts3D; FRAMES = FRAMES_3D;
% Field to plot
% FIELD = dens_e;       NAME = 'ne';   FIELDNAME = 'n_e';
% FIELD = dens_i;       NAME = 'ni';   FIELDNAME = 'n_i';
% FIELD = dens_e-Z_n_e; NAME = 'ne_NZ';FIELDNAME = 'n_e^{NZ}';
% FIELD = dens_i-Z_n_i; NAME = 'ni_NZ';FIELDNAME = 'n_i^{NZ}';
% FIELD = temp_e;       NAME = 'Te';   FIELDNAME = 'n_i';
% FIELD = temp_i;       NAME = 'Ti';   FIELDNAME = 'n_i';
% FIELD = temp_e-Z_T_e; NAME = 'Te_NZ';FIELDNAME = 'T_e^{NZ}';
% FIELD = temp_i-Z_T_i; NAME = 'Ti_NZ';FIELDNAME = 'T_i^{NZ}';
% FIELD = ne00;         NAME = 'ne00'; FIELDNAME = 'n_e^{00}';
% FIELD = ni00;   NAME = 'ni00'; FIELDNAME = 'n_i^{00}';
% FIELD = phi;    NAME = 'phi'; FIELDNAME = '\phi';
% FIELD = Gamma_x;  NAME = 'Gamma_x'; FIELDNAME = '\Gamma_x';

% Sliced
% plt = @(x) real(x(ix, :, :,:)); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; 
% plt = @(x) real(x( :,iy, :,:)); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; 
% plt = @(x) real(x( :, :,iz,:)); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; 

% Averaged
% plt = @(x) mean(x,1); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; 
% plt = @(x) mean(x,2); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; 
plt = @(x) mean(x,3); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; 

FIELD = squeeze(plt(FIELD));

% Naming
GIFNAME   = [NAME,sprintf('_%.2d',JOBNUM),'_',PARAMS];

% Create movie (gif or mp4)
create_gif
% create_mov
end

if 0
%% Photomaton : real space

% Chose the field to plot
% FIELD = ni00;          FNAME = 'ni00';    FIELDLTX = 'n_i^{00}';
% FIELD = ne00;          FNAME = 'ne00';    FIELDLTX = 'n_e^{00}'
% FIELD = dens_i;        FNAME = 'ni';      FIELDLTX = 'n_i';
% FIELD = dens_e;        FNAME = 'ne';      FIELDLTX = 'n_e';
FIELD = dens_e-Z_n_e;   FNAME = 'ne_NZ';   FIELDLTX = 'n_e^{NZ}';
% FIELD = dens_i-Z_n_i;   FNAME = 'ni_NZ';   FIELDLTX = 'n_i^{NZ}';
% FIELD = temp_i;        FNAME = 'Ti';      FIELDLTX = 'T_i';
% FIELD = temp_e;        FNAME = 'Te';      FIELDLTX = 'T_e';
% FIELD = phi;           FNAME = 'phi';     FIELDLTX = '\phi';
% FIELD = Z_phi-phi;     FNAME = 'phi_NZ';  FIELDLTX = '\phi^{NZ}';
% FIELD = Gamma_x;       FNAME = 'Gamma_x'; FIELDLTX = '\Gamma_x';
% FIELD = dens_e-Z_n_e-(Z_phi-phi);       FNAME = 'Non_adiab_part'; FIELDLTX = 'n_e^{NZ}-\phi^{NZ}';

% Chose when to plot it
tf = 500:500:2500;

% Sliced
ix = 1; iy = 1; iz = 1;
% plt = @(x,it) real(x(ix, :, :,it)); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; FIELDLTX = [FIELDLTX,'(x=',num2str(round(x(ix))),')']
% plt = @(x,it) real(x( :,iy, :,it)); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; FIELDLTX = [FIELDLTX,'(y=',num2str(round(y(iy))),')']
plt = @(x,it) real(x( :, :,iz,it)); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; FIELDLTX = [FIELDLTX,'(z=',num2str(round(z(iz)/pi)),'\pi)'] 

% Averaged
% plt = @(x,it) mean(x(:,:,:,it),1); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; FIELDLTX = ['\langle ',FIELDLTX,'\rangle_x']
% plt = @(x,it) mean(x(:,:,:,it),2); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; FIELDLTX = ['\langle ',FIELDLTX,'\rangle_y']
% plt = @(x,it) mean(x(:,:,:,it),3); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; FIELDLTX = ['\langle ',FIELDLTX,'\rangle_z'] 


%
TNAME = [];
fig = figure; FIGNAME = [FNAME,TNAME,'_snaps','_',PARAMS]; set(gcf, 'Position',  [100, 100, 1500, 350])
plt_2 = @(x) x./max(max(x));
    for i_ = 1:numel(tf)
    [~,it] = min(abs(Ts3D-tf(i_))); TNAME = [TNAME,'_',num2str(Ts3D(it))];
    subplot(1,numel(tf),i_)
        DATA = plt_2(squeeze(plt(FIELD,it)));
        pclr = pcolor((X),(Y),DATA); set(pclr, 'edgecolor','none');pbaspect([1 1 1])
        colormap(bluewhitered); caxis([-1,1]);
        xlabel(latexize(XNAME)); ylabel(latexize(YNAME));set(gca,'ytick',[]); 
        title(sprintf('$t c_s/R=%.0f$',Ts3D(it)));
    end
    legend(latexize(FIELDLTX));
save_figure
end

if 0
%% Photomaton : k space

% Chose the field to plot
% FIELD = Ni00;   FNAME = 'Ni00'; FIELDLTX = 'N_i^{00}';
FIELD = Ne00;   FNAME = 'Ne00'; FIELDLTX = 'N_e^{00}'
% FIELD = PHI; FNAME = 'PHI'; FIELDLTX = '\tilde\phi';
% FIELD_ = fft2(Gamma_x); FIELD = FIELD_(1:Nx/2+1,:,:,:); FNAME = 'Gamma_x'; FIELDLTX = '\tilde\Gamma_x';
% FIELD_ = fft2(dens_e); FIELD = FIELD_(1:Nx/2+1,:,:,:); FNAME = 'FFT_Dens_e'; FIELDLTX = '\tilde n_e';

% Chose when to plot it
tf = 500:500:2500;

% Sliced
ix = 1; iy = 1; iz = 1;
% plt = @(x,it) abs(x(ix, :, :,it)); X = KY_YZ; Y = KZ_YZ; XNAME = 'k_y'; YNAME = 'z'; FIELDLTX = [FIELDLTX,'(k_x=',num2str(round(kx(ix))),')'];
% plt = @(x,it) abs(x( :,iy, :,it)); X = KX_XZ; Y = KZ_XZ; XNAME = 'k_x'; YNAME = 'z'; FIELDLTX = [FIELDLTX,'(k_y=',num2str(round(ky(iy))),')'];
plt = @(x,it) abs(x( :, :,iz,it)); X = KX_XY; Y = KY_XY; XNAME = 'k_x'; YNAME = 'k_y'; FIELDLTX = [FIELDLTX,'(z=',num2str((z(iz)/pi)),'\pi)']; 

%
TNAME = [];
fig = figure; FIGNAME = [FNAME,TNAME,'_snaps','_',PARAMS]; set(gcf, 'Position',  [100, 100, 300*numel(tf), 500])
plt_2 = @(x) (fftshift(x,2));
    for i_ = 1:numel(tf)
    [~,it] = min(abs(Ts3D-tf(i_))); TNAME = [TNAME,'_',num2str(Ts3D(it))];
    subplot(1,numel(tf),i_)
        DATA = plt_2(squeeze(plt(FIELD,it)));
        pclr = pcolor(fftshift(X,2),fftshift(Y,2),DATA); set(pclr, 'edgecolor','none');pbaspect([0.5 1 1])
%         colormap(bluewhitered); caxis([-1,1]);
        xlabel(latexize(XNAME)); ylabel(latexize(YNAME));
        if(i_ > 1); set(gca,'ytick',[]); end;
        title(sprintf('$t c_s/R=%.0f$',Ts3D(it)));
    end
    legend(latexize(FIELDLTX));
save_figure
end

%%
% ZF_fourier_analysis