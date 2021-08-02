addpath(genpath('../matlab')) % ... add
addpath(genpath('../matlab/plots')) % ... add
%% Directory of the simulation
if 1% Local results
outfile ='HD_study/150x75_L_100_P_4_J_2_eta_0.6_nu_1e-02_DGGK_mu_3e-03';
    BASIC.RESDIR      = ['../results/',outfile,'/'];
    BASIC.MISCDIR     = ['/misc/HeLaZ_outputs/results/',outfile,'/'];
    CMD = ['cp ', BASIC.RESDIR,'outputs* ',BASIC.MISCDIR]; disp(CMD);
    system(CMD);
end
if 0% Marconi results
outfile ='';
outfile ='';
outfile ='';
outfile ='';
outfile ='';
outfile ='';
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

if 1
%% Time evolutions and growth rate
plot_time_evolution_and_gr
end

if 1
%% Space time diagramm (fig 11 Ivanov 2020)
TAVG = 5000; % Averaging time duration
plot_radial_transport_and_shear
end

if 0
%% Space time diagramms
cmax = 0.01 % max of the colorbar for transport
tstart = 0; tend = Ts3D(end); % time window
plot_space_time_diagrams
end

if 0
%% Averaged shear and Reynold stress profiles
trange = its2D:ite2D;
plot_shear_and_reynold_stress
end

if 0
%% |phi_k|^2 spectra (Kobayashi 2015 fig 3)
tstart = 0.8*Ts3D(end); tend = Ts3D(end); % Time window
plot_kperp_spectrum
end

if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
t0    =00; iz = 1; ix = 1; iy = 1;
skip_ = 1; DELAY = 1e-2*skip_;
[~, it03D] = min(abs(Ts3D-t0)); FRAMES_3D = it03D:skip_:numel(Ts3D);
[~, it05D] = min(abs(Ts5D-t0)); FRAMES_5D = it05D:skip_:numel(Ts5D);
INTERP = 0; T = Ts3D; FRAMES = FRAMES_3D;
% Field to plot
FIELD = dens_e; NAME = 'ne';   FIELDNAME = 'n_e';
% FIELD = dens_i; NAME = 'ni';   FIELDNAME = 'n_i';
% FIELD = temp_e; NAME = 'Te';   FIELDNAME = 'n_i';
% FIELD = temp_i; NAME = 'Ti';   FIELDNAME = 'n_i';
% FIELD = ne00;   NAME = 'ne00'; FIELDNAME = 'n_e^{00}';
% FIELD = ni00;   NAME = 'ni00'; FIELDNAME = 'n_i^{00}';
% Slice
% plt = @(x) real(x(ix, :, :,:)); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; 
% plt = @(x) real(x( :,iy, :,:)); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; 
plt = @(x) real(x( :, :,iz,:)); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; 

% Averaged
% plt = @(x) mean(x,1); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; 
% plt = @(x) mean(x,2); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; 
% plt = @(x) mean(x,3); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; 

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
% FIELD = ni00;   FNAME = 'ni00'; FIELDLTX = 'n_i^{00}';
% FIELD = ne00;   FNAME = 'ne00'; FIELDLTX = 'n_e^{00}'
% FIELD = dens_i; FNAME = 'ni';   FIELDLTX = 'n_i';
% FIELD = dens_e; FNAME = 'ne';   FIELDLTX = 'n_e';
% FIELD = temp_i; FNAME = 'Ti';   FIELDLTX = 'T_i';
% FIELD = temp_e; FNAME = 'Te';   FIELDLTX = 'T_e';
FIELD = phi; FNAME = 'phi'; FIELDLTX = '\phi';

% Chose when to plot it
tf = [0 1 2 3];

% Slice
ix = 1; iy = 1; iz = 1;
% plt = @(x,it) real(x(ix, :, :,it)); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; FIELDLTX = [FIELDLTX,'(x=',num2str(round(x(ix))),')']
% plt = @(x,it) real(x( :,iy, :,it)); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; FIELDLTX = [FIELDLTX,'(y=',num2str(round(y(iy))),')']
% plt = @(x,it) real(x( :, :,iz,it)); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; FIELDLTX = [FIELDLTX,'(x=',num2str(round(z(iz))),')'] 

% Averaged
% plt = @(x,it) mean(x(:,:,:,it),1); X = Y_YZ; Y = Z_YZ; XNAME = 'y'; YNAME = 'z'; FIELDLTX = ['\langle',FIELDLTX,'\rangle_x']
% plt = @(x,it) mean(x(:,:,:,it),2); X = X_XZ; Y = Z_XZ; XNAME = 'x'; YNAME = 'z'; FIELDLTX = ['\langle',FIELDLTX,'\rangle_y']
plt = @(x,it) mean(x(:,:,:,it),3); X = X_XY; Y = Y_XY; XNAME = 'x'; YNAME = 'y'; FIELDLTX = ['\langle',FIELDLTX,'\rangle_z'] 


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

%%
% ZF_fourier_analysis