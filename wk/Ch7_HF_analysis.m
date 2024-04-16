PARTITION = '/misc/gyacomo23_outputs/triangularity_paper/';
switch 4
case 1
    SIM_SET_NAME = 'Multi-scale';
    E_FLUX       = 1;
    FOLDER       = 'DIIID_fullphys_rho95_geom_scan/multi_scale/multi_scale_3x2x768x192x24';
case 2
    SIM_SET_NAME = 'KEM';
    E_FLUX       = 1;
    FOLDER       = 'ion_scale/5x2x256x64x32';
case 3
    SIM_SET_NAME = 'AEM';
    E_FLUX       = 0;
    FOLDER       = 'adiabatic_electrons/5x2x256x64x32';
case 4
    SIM_SET_NAME = 'RFM';
    E_FLUX       = 0;
    FOLDER       = 'hot_electrons/256x64x32';
end

GEOM = {'NT','0T','PT'};
COLO = {'b','k','r'};
figure
for i = 1:3
    OPTION = GEOM{i};
    datadir = [PARTITION,FOLDER,'/',OPTION];
    out = read_flux_out_XX(datadir,0);
    
    [ts, is] = sort(out.t);
    Pxis = out.Pxi(is);
    Qxis = out.Qxi(is);

    if E_FLUX
        Pxes = out.Pxe(is);
        Qxes = out.Qxe(is);
        subplot(211)
        plot(ts,Qxes,COLO{i}); hold on;
        subplot(212)
        title(SIM_SET_NAME)
    end
    plot(ts,Qxis,COLO{i}); hold on;
end
title(SIM_SET_NAME)

