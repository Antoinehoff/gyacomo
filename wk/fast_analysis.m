% Directory of the code "mypathtogyacomo/gyacomo/"
% Partition of the computer where the data have to be searched
PARTITION  = '/misc/gyacomo_outputs/';

%% Scan kT
resdirs = {...
'paper_2_nonlinear/kT_scan_nu_1e-3/5x3x128x64x24_dp', ...
'paper_2_nonlinear/kT_scan_nu_1e-3/5x3x192x96x32_dp', ...
'paper_2_nonlinear/kT_scan_nu_1e-3/7x4x128x64x24_dp', ...
'paper_2_nonlinear/kT_scan_nu_1e-3/7x4x192x96x32_dp', ...
};

%% Scan nu, kT = 5.3
% resdirs = {...
% 'paper_2_nonlinear/nu_scan_kT_5.3/FCGK_5x3x128x64x24_dp', ...
% 'paper_2_nonlinear/nu_scan_kT_5.3/DGGK_7x4x128x64x24_dp', ...
% 'paper_2_nonlinear/nu_scan_kT_5.3/SGGK_7x4x128x64x24_dp', ...
% };


%%
figure
hold on
for i = 1:numel(resdirs)
    J0 = 00; J1 = 10;

    % Load basic info (grids and time traces)
    DATADIR = [PARTITION,resdirs{i},'/'];
    data    = {};
    data    = compile_results_low_mem(data,DATADIR,J0,J1);
    
    % plot heat flux
    subplot(1,2,1)
    hold on
    plot(data.Ts0D,data.HFLUX_X,'DisplayName',data.paramshort);
    
    % statistical transport averaging
    Gavg =[]; Gstd = [];
    Qavg =[]; Qstd = [];
    for i_ = 1:2:numel(data.TJOB_SE) 
    disp([num2str(data.TJOB_SE(i_)),' ',num2str(data.TJOB_SE(i_+1))])
    disp([num2str(data.NU_EVOL(i_)),' ',num2str(data.NU_EVOL(i_+1))])
    options.T = [data.TJOB_SE(i_)*1.2 data.TJOB_SE(i_+1)];
    options.NPLOTS = 0;
    [fig, res] = statistical_transport_averaging(data,options);
    Gavg = [Gavg res.Gx_avg]; Gstd = [Gstd res.Gx_std];
    Qavg = [Qavg res.Qx_avg]; Qstd = [Qstd res.Qx_std];
    end
    subplot(1,2,2)
    hold on
    errorbar(data.K_T_EVOL(2:2:end),Qavg,Qstd,'--s','DisplayName',data.paramshort);xlabel('$\kappa_T$'); 
%     errorbar(data.NU_EVOL(2:2:end),Qavg,Qstd,'--s','DisplayName',data.paramshort);xlabel('$\nu$'); 
    ylabel('$Q_x^\infty$');
end

if 0
%% Plot transport and phi radial profile
[data.PHI, data.Ts3D] = compile_results_3D(DATADIR,J0,J1,'phi');

options.TAVG_0   = 100;
options.TAVG_1   = 1000;
options.NCUT     = 5;              % Number of cuts for averaging and error estimation
options.NMVA     = 1;              % Moving average for time traces
% options.ST_FIELD = '\Gamma_x';   % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.ST_FIELD = '\phi';          % chose your field to plot in spacetime diag (e.g \phi,v_x,G_x)
options.INTERP   = 0;
options.RESOLUTION = 256;
fig = plot_radial_transport_and_spacetime(data,options);
end
if 0
%% MOVIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
options.INTERP    = 0;
options.POLARPLOT = 0;
options.NAME      = '\phi';
% options.NAME      = '\omega_z';
% options.NAME     = 'N_i^{00}';
% options.NAME      = 's_{Ey}';
% options.NAME      = 'n_i^{NZ}';
% options.NAME      = 'Q_x';
% options.NAME      = 'n_i';
% options.NAME      = 'n_i-n_e';
options.PLAN      = 'xz';
% options.NAME      = 'f_i';
% options.PLAN      = 'sx';
options.COMP      = 'avg';
% options.TIME      = data.Ts5D(end-30:end);
% options.TIME      =  data.Ts3D;
options.TIME      = [0:1500];
data.EPS          = 0.1;
data.a = data.EPS * 2000;
options.RESOLUTION = 256;
create_film(data,options,'.gif')
end